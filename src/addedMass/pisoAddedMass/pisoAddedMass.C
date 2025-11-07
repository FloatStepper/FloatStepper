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

#include "pisoAddedMass.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "slipVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(pisoAddedMass, 0);
    addToRunTimeSelectionTable
    (
        addedMass,
        pisoAddedMass,
        dict
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pisoAddedMass::pisoAddedMass
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh)
{}

Foam::SquareMatrix<Foam::scalar> Foam::pisoAddedMass::computeAddedMass
(
    FixedList<bool,6> enabledDoF,
    const vector& CoR,
    const tensor& Q,
    const wordRes patches
)
{

    const scalar deltaT = mesh_.time().deltaTValue();

    // Setting size of virtual linear and angular accelerations
    const scalar dv = 1;
    const scalar domega = 1;
    const scalar linAcc = dv/deltaT;
    const scalar rotAcc = domega/deltaT;
    Info << "Probe accelerations: linAcc = " << linAcc << ", rotAcc = " 
        << rotAcc << endl;

    // Make a loop around the section below to populate added mass matrix
    typedef SquareMatrix<scalar> SMatrix;
    SMatrix Madd({6, 6}, 0); // Added mass tensor in lab frame

    forAll(enabledDoF, DoFi)
    {
        if (enabledDoF[DoFi])
        {
            Info << "------------Virtual displacement for DoF " << DoFi
                << "-----------" << endl;
            vector dvVec(0, 0, 0);
            vector dwVec(0, 0, 0);
            if (DoFi < 3) dvVec[DoFi] = dv;
            else dwVec[DoFi-3] = domega;
            vector F = Zero;
            vector tau = Zero;
            calcPressureForceAndTorque(mesh_, dvVec, dwVec, CoR, patches, F, tau);

            // Recording fluid force associated with small force change motion
            Info << "DoF = " << DoFi << ", dvVec = " << dvVec << ", dwVec = "
                << dwVec << endl;
            Info << "F = " << F << ", tau = " << tau << endl;
            Madd.subMatrix(0, DoFi, 3, 1) = - F/linAcc;
            Madd.subMatrix(3, DoFi, 3, 1) = - tau/rotAcc;
        }
    }

    return Madd;
}


void Foam::pisoAddedMass::calcPressureForceAndTorque
(
    const fvMesh& mesh,
    const vector& v,
    const vector& omega,
    const point& CoR,
    const wordRes& patches,
    vector& force,
    vector& torque
)
{
    // Reading fields that we must copy to take added mass PISO step
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const surfaceScalarField& phi = mesh.lookupObject<surfaceScalarField>("phi");
    const volScalarField& p_rgh = mesh.lookupObject<volScalarField>("p_rgh");
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");

    const pimpleControl& pimple = mesh.lookupObject<pimpleControl>("solutionControl");

    const auto& bMesh = mesh.boundary();

    // Copy of U to get all its boundary conditions identical
    volVectorField U2("U2", U);
    U2.oldTime() = (0*U2);

    // Overwriting floating object boundary conditions for U2 to correspond to
    // rigid body slip conditions.
    forAll(patches, patchi)
    {
        const label patchId
            = mesh.boundaryMesh().findPatchID(patches[patchi], false);

        U2.boundaryFieldRef().set
        (
            patchId,
            fvPatchField<vector>::New
            (
                "slipVelocity",
                bMesh[patchId],
                U2
            )
        );

        refCast<slipVelocityFvPatchVectorField>
        (
            U2.boundaryFieldRef()[patchId]
        ).updateCoeffs(CoR, v, omega);
    }

    surfaceScalarField phi2(phi);
    // Note: Important that p_rgh2 is given a unique name
    // Otherwise force and torque will be calculated incorrectly.
    volScalarField p_rgh2("p_rgh2", p_rgh);
    p_rgh2 = dimensionedScalar("zero", p_rgh.dimensions(), Zero);
    p_rgh2.oldTime() = (0*p_rgh2);
    mesh.setFluxRequired(p_rgh2.name());

    fvVectorMatrix U2Eqn(fvm::ddt(rho, U2));

    // --- Pressure corrector loop
    for (int m = 1; m <= pimple.nCorrPISO(); m++)
    {

        volScalarField rAU2(1.0/U2Eqn.A());
        surfaceScalarField rAU2f("rAU2f", fvc::interpolate(rAU2));

        volVectorField HbyA2(constrainHbyA(rAU2*U2Eqn.H(), U2, p_rgh2));
        surfaceScalarField phiHbyA2("phiHbyA2", fvc::flux(HbyA2));

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh2, U2, phiHbyA2, rAU2f);
        // Question: Why not use rho in this expression?: 
        // constrainPressure(p_rgh2, rho, U2, phiHbyA2, rAU2f);

        label pRefCell = 0;
        scalar pRefValue = 0.0;
        setRefCell(p, p_rgh2, pimple.dict(), pRefCell, pRefValue);

        for (label n = 0; n < pimple.nNonOrthCorr() + 1 ; n++)
        {
            fvScalarMatrix p_rgh2Eqn
            (
                fvm::laplacian(rAU2f, p_rgh2) == fvc::div(phiHbyA2)
            );

            p_rgh2Eqn.setReference(pRefCell, getRefCellValue(p_rgh2, pRefCell));

            p_rgh2Eqn.solve
            (
                p_rgh2.mesh().solver(p_rgh.select(pimple.finalInnerIter()))
            );

            if (n == pimple.nNonOrthCorr())
            {
                phi2 = phiHbyA2 - p_rgh2Eqn.flux();

                U2 = HbyA2 - rAU2*fvc::reconstruct(p_rgh2Eqn.flux()/rAU2f);
                U2.correctBoundaryConditions();
            }
        }
    }

    dictionary forcesDict;

    forcesDict.add("type", functionObjects::forces::typeName);
    forcesDict.add("patches", patches);
    forcesDict.add("rhoInf", 1.0);
    forcesDict.add("rho", rho.name());
    forcesDict.add("p", p_rgh2.name());
    forcesDict.add("U", U2.name());
    forcesDict.add("CofR", CoR);

    functionObjects::forceContributions f("forces", rho.time(), forcesDict);

    f.calcForcesMoments();
    force = f.pressureForce();
    torque = f.pressureMoment();
}


// ************************************************************************* //
