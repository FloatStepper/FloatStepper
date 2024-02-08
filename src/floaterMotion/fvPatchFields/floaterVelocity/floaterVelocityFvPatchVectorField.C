/*---------------------------------------------------------------------------*\
    Module Name:       FloatStepper
    Description:       OpenFOAM extension module for fluid-rigid body coupling
    License:           GNU General Public License (GPL) version 3
    Copyright:         2023 Johan Roenby, STROMNING APS
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of FloatStepper.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FloatStepper.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "floaterVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
//#include "fvcMeshPhi.H"
#include "symmTransformField.H"

#include "motionSolver.H"
#include "floaterMotionSolver.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterVelocityFvPatchVectorField::
floaterVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    bodyName_(),
    slip_(false)
{}


Foam::floaterVelocityFvPatchVectorField::
floaterVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    bodyName_(),
    slip_(dict.lookupOrDefault<bool>("slip", false))
{

    // Setting bodyName_ by looking up the motionSolver of type floaterMotion
    // in dynamicMeshDict.solvers that has this patch in its patches list.
    IOobject ioDict
    (
        "dynamicMeshDict",
        db().time().constant(),
        patch().boundaryMesh().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    IOdictionary dynMeshDict(ioDict);

    if (dynMeshDict.found("solvers"))
    {
        const dictionary& solvertDict = dynMeshDict.subDict("solvers");

        for (const entry& dEntry : solvertDict)
        {
            if (dEntry.isDict() && dEntry.dict().found("motionSolver"))
            {
                const word ms(dEntry.dict().lookup("motionSolver"));
                if ( ms == "floaterMotion" )
                {
                    const dictionary& RBMDict = dEntry.dict().subDict(ms+"Coeffs");
                    const wordRes patchNames(RBMDict.get<wordRes>("patches"));
                    forAll(patchNames, ni)
                    {
                        if (patchNames[ni] == patch().name())
                        {
                            bodyName_ = dEntry.dict().dictName();
                        }
                    }
                }
            }
        }
    }
    else
    {
        bodyName_ = dynMeshDict.dictName();
        WarningInFunction
        << "Subdict named solvers not found in " << bodyName_ << "."
        << " . Setting bodyName to " << bodyName_ << " on patch "
        << patch().name() << endl;
    }

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        // Evaluate the wall velocity
        updateCoeffs();
    }

}


Foam::floaterVelocityFvPatchVectorField::
floaterVelocityFvPatchVectorField
(
    const floaterVelocityFvPatchVectorField& fvfpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(fvfpvf, p, iF, mapper),
    bodyName_(fvfpvf.bodyName_),
    slip_(fvfpvf.slip_)
{}


Foam::floaterVelocityFvPatchVectorField::
floaterVelocityFvPatchVectorField
(
    const floaterVelocityFvPatchVectorField& fvfpvf
)
:
    fixedValueFvPatchField<vector>(fvfpvf),
    bodyName_(fvfpvf.bodyName_),
    slip_(fvfpvf.slip_)
{}


Foam::floaterVelocityFvPatchVectorField::
floaterVelocityFvPatchVectorField
(
    const floaterVelocityFvPatchVectorField& fvfpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(fvfpvf, iF),
    bodyName_(fvfpvf.bodyName_),
    slip_(fvfpvf.slip_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::floaterVelocityFvPatchVectorField::updateCoeffs()
{

    if (updated())
    {
        return;
    }

    vector v(Zero);
    vector omega(Zero);
    vector CoR(Zero);

    const fvMesh& mesh = this->internalField().mesh();

    if (mesh.foundObject<floaterMotionSolver>(bodyName_))
    {
        const floaterMotionSolver& bodySolver =
            mesh.lookupObject<floaterMotionSolver>(bodyName_);

        const floaterMotionState& bodyState(bodySolver.motion().state());

        v = bodyState.v();
        omega = bodyState.omega();
        CoR = bodyState.centreOfRotation();
    }
    else
    {
        WarningInFunction
        << "floaterMotionSolver named " << bodyName_ 
        << " not found in objectRegistry."
        << "Velocity set to zero on patch " << patch().name() << endl;
    }

    const vectorField Up(v + (omega^(patch().Cf() - CoR)));

    if (slip_)
    {
        // Note: Using reference here causes crash
        const vectorField nHat(this->patch().nf());

        vectorField::operator=
        (
            (nHat & Up)*nHat
                + transform(I - sqr(nHat), this->patchInternalField())
        );
    }
    else
    {
        vectorField::operator=(Up);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::floaterVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("bodyName", bodyName_);
    os.writeEntry("slip", slip_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        floaterVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
