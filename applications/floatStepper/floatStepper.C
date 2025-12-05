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
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 DHI
    Copyright (C) 2017 OpenCFD Ltd.
    Copyright (C) 2018 Johan Roenby
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2020 OpenCFD Ltd.
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

Application
    floatStepper

Description
    Solver for coupled motion of a rigid floating body and a surrounding fluid
    consisting of two incompressible, isothermal immiscible fluids using a
    geometric VOF (isoAdvector) phase-fraction based interface capturing
    approach with mesh morphing to accommodate the body motion.

    Johan Roenby, STROMNING 2023

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "floaterMeshMotionSolver.H"
#include "myAdjustPhi.H"
#include "myCorrectPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for rigid body floating in two incompressible, isothermal "
        "immiscible fluids. Uses the isoAdvector geometric VOF method for "
        "interface capturing with mesh morphing to accommodate the body motion."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        const scalar deltaT = runTime.deltaTValue();

        List<floaterMotionState> oldBodyStates(bodySolvers.size());

        forAll(bodySolvers, bsi)
        {
            floaterMeshMotionSolver& bodySolver = bodySolvers[bsi];

            if (onlyMeshMotion)
            {
                Info << "Only moving mesh - no body state update." << endl;
                vector a = bodySolver.motion().state().a();
                vector alpha = bodySolver.motion().state().domegadt();
                scalarField dvwdt(6, 0);
                dvwdt[0] = a[0], dvwdt[1] = a[1], dvwdt[2] = a[2];
                dvwdt[3] = alpha[0], dvwdt[4] = alpha[1], dvwdt[5] = alpha[2];
                bodySolver.motion().updateFloaterState(dvwdt, deltaT);
                //Info << "Body state after first time step:" << endl;
                //bodySolver.motion().status();
            }
            else
            {

                // Copying point fields for resetting.
                //pointField oldPoints(mesh.points());
                //pointField oldOldPoints(mesh.oldPoints());

                // Saving previous time step body state
                //const floaterMotionState oldBodyState(bodySolver.motion().state());
                oldBodyStates[bsi] = bodySolver.motion().state();
                // Calculating the 6 columns of the added mass tensor
                const labelList DoFs = bodySolver.DoFs();

                if (bodySolver.motion().MaddUpdateTime())
                {
                    mesh.moving(false);
                    bodySolver.motion().calcAddedMass(mesh, DoFs);
                    mesh.moving(true);
                }

                // Taking zero acceleration time step to get F0 and tau0
                const scalarField noAcceleration(6,0);
                bodySolver.motion().setAcceleration(Zero);
                bodySolver.motion().setAngularAcceleration(Zero);
                bodySolver.motion().updateFloaterState(noAcceleration, deltaT);
                //Info << "Body state after zero acceleration time step:" << endl;
                //bodySolver.motion().status();
            }
        }

        // Now all body states are evolved with zero acceleration

        // Take zero acceleration time step
        vectorList F0List(bodySolvers.size(), Zero);
        vectorList tau0List(bodySolvers.size(), Zero);
        {
            // Update mesh points in accordance with 0-acceleration body motion
            #include "updateMesh0.H"

            // Calculate fluid motion response to 0-acceleration body motion
            #include "updateFluid0.H" //alphaEqn, UEqn and Piso loop

            // Recording fluid force associated with 0-acceleration motion
            forAll(bodySolvers, bsi)
            {
                floaterMeshMotionSolver& bodySolver = bodySolvers[bsi];
                vector F0(Zero), tau0(Zero);
                bodySolver.motion().calcForceAndTorque(rho0, p0, U0, F0, tau0);
                F0List[bsi] = F0; tau0List[bsi] = tau0;
                Info << "F0 = " << F0 << ", tau0 = " << tau0 << endl;
            }

            // Reset fluid and body states
            #include "reset.H"
        }

        // Calculate acceleration and update state for each floater
        if (!onlyMeshMotion)
        {
            forAll(bodySolvers, bsi)
            {
                floaterMeshMotionSolver& bodySolver = bodySolvers[bsi];
                const labelList DoFs = bodySolver.DoFs();

                // Update body state based on added mass, F0 and tau0
                scalarField dvwdt =
                    bodySolver.motion().calcAcceleration(F0List[bsi], tau0List[bsi], DoFs);
                // Updating body state to new time
                bodySolver.motion().updateFloaterState(dvwdt, deltaT);
                Info << "Body state after real time step:" << endl;
                bodySolver.motion().status();
                Info << "Centre of mass      : "
                    << bodySolver.motion().centreOfMass() << endl;
            }
        }

        // Updating mesh and fluid state corresponding to body motion
        Info << "Taking actual time step for fluid..." << endl;
        #include "updateMesh.H"
        #include "updateFluid.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
