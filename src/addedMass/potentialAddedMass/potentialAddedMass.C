/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-----------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2023 Johan Roenby
-----------------------------------------------------------------------------
License
    This file is part of FloatStepper.

    FloatStepper is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FloatStepper is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FloatStepper.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "potentialAddedMass.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "rigidBodyVelocityPotentialFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(potentialAddedMass, 0);
    addToRunTimeSelectionTable
    (
        addedMass,
        potentialAddedMass,
        dict
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::potentialAddedMass::potentialAddedMass
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Pot_
    (
        IOobject
        (
            "Pot",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Pot", dimless, 0),
        "zeroGradient"
    )
{}

Foam::SquareMatrix<Foam::scalar> Foam::potentialAddedMass::computeAddedMass
(
    FixedList<bool,6> enabledDoF,
    const vector& CoR,
    const tensor& Q,
    const wordRes patches
)
{
    const fvMesh& mesh = Pot_.mesh();

    const auto& bMesh = mesh.boundary();

    Info << "enabledDoF " << enabledDoF << endl;

    const auto& Potf = Pot_.boundaryField();
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    const auto& rhof = rho.boundaryField();

    typedef SquareMatrix<scalar> SMatrix;
    SMatrix Madd({6, 6}, 0);

    forAll(enabledDoF, DoF)
    {
        bool validDir = enabledDoF[DoF];
        if (validDir)
        {
            Pot_ *= 0.0; // reset Pot

            vector v = Zero;
            vector omega = Zero;
            if (DoF < 3) v[DoF] = 1;
            else omega[DoF-3] = 1;

            forAll(patches, patchi)
            {
                // fatal error if patch does not excist
                const label patchId
                    = mesh.boundaryMesh().findPatchID(patches[patchi], false);

                Pot_.boundaryFieldRef().set
                (
                    patchId,
                    fvPatchField<scalar>::New
                    (
                        "rigidBodyVelocityPotential",
                        bMesh[patchId],
                        Pot_
                    )
                );

                const vectorField nf(bMesh[patchId].nf());
                const vectorField Cf(bMesh[patchId].Cf());

                refCast<rigidBodyVelocityPotentialFvPatchScalarField>
                (
                    Pot_.boundaryFieldRef()[patchId]
                ).updateSnGrad((nf & (v + (omega ^ (Cf - CoR)))));

            }

            fvScalarMatrix potEqn(fvm::laplacian(rho, Pot_));

            potEqn.solve();

            forAll(patches, patchi)
            {
                const label patchId
                    = mesh.boundaryMesh().findPatchID(patches[patchi], false);

                const vectorField Cf(bMesh[patchId].Cf());
                const vectorField Sf(bMesh[patchId].Sf());

                SMatrix Maddp({6, 6}, 0);

                Maddp.subMatrix(0, DoF, 3, 1) 
                    = gSum(rhof[patchId]*Potf[patchId]*Sf);
                Maddp.subMatrix(3, DoF, 3, 1) 
                    = gSum(rhof[patchId]*Potf[patchId]*((Cf - CoR) ^ Sf));

                Madd += Maddp;
            }
        }
    }

    return Madd;
}



// ************************************************************************* //
