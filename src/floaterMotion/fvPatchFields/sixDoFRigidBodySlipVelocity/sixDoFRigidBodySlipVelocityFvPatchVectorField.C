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

#include "sixDoFRigidBodySlipVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "symmTransformField.H"

#include "motionSolver.H"
#include "sixDoFRigidBodyMotionSolver.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodySlipVelocityFvPatchVectorField::
sixDoFRigidBodySlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{}


Foam::sixDoFRigidBodySlipVelocityFvPatchVectorField::
sixDoFRigidBodySlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false)
{
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


Foam::sixDoFRigidBodySlipVelocityFvPatchVectorField::
sixDoFRigidBodySlipVelocityFvPatchVectorField
(
    const sixDoFRigidBodySlipVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::sixDoFRigidBodySlipVelocityFvPatchVectorField::
sixDoFRigidBodySlipVelocityFvPatchVectorField
(
    const sixDoFRigidBodySlipVelocityFvPatchVectorField& rwvpvf
)
:
    fixedValueFvPatchField<vector>(rwvpvf)
{}


Foam::sixDoFRigidBodySlipVelocityFvPatchVectorField::
sixDoFRigidBodySlipVelocityFvPatchVectorField
(
    const sixDoFRigidBodySlipVelocityFvPatchVectorField& rwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(rwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodySlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();

    if (mesh.foundObject<sixDoFRigidBodyMotionSolver>("dynamicMeshDict"))
    {

        sixDoFRigidBodyMotionSolver& bodySolver =
            const_cast<sixDoFRigidBodyMotionSolver&>
            (
                mesh.lookupObject<sixDoFRigidBodyMotionSolver>("dynamicMeshDict")
            );
        const sixDoFRigidBodyMotionState& bodyState(bodySolver.motion().state());

        vector v = bodyState.v();
        vector omega = bodySolver.motion().omega();
        vector CoR = bodyState.centreOfRotation();

        const vectorField nHat(this->patch().nf());

        const vectorField Up
        (
            (nHat & (v + (omega^(patch().Cf() - CoR))))*nHat
                + transform(I - sqr(nHat), this->patchInternalField())
        );

        vectorField::operator=(Up);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::sixDoFRigidBodySlipVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        sixDoFRigidBodySlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
