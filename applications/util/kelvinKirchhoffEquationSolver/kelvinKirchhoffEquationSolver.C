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
    kirchhoffKelvinEquationSolver

Description
    Solves the 6-DoF Kelvin-Kirchhoff equations of motion for a rigid body in an
    ideal fluid (also sometimes called the Kirchhoff-Kelvin equations).

    Johan Roenby, STROMNING 2025

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "subCycle.H"
#include "fvOptions.H"
#include "floaterMotion.H"
#include "ODESystem.H"
#include "ODESolver.H"

using namespace Foam;
using namespace Foam::constant::mathematical;
using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "kelvinKirchhoffEquations.H"

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
    #include "createFields.H"

    // Read floaters file
    IOdictionary floaterStateDict
    (
        IOobject
        (
            "floaters",
            mesh.time().timeName(),
            "uniform",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    // Get the list of all bodies
    wordList bodies = floaterStateDict.toc();
    word bodyName;
    // Pick first body in file
    forAll(bodies, i)
    {
        const word& bodyi = bodies[i];
        if (floaterStateDict.isDict(bodyi))
        {
            bodyName = bodyi;
            break;
        }
    }
    if (bodies.size() > 1)
    {
        Info << "More than one body in floaters file" << endl;
        Info << "Choosing the first one named: " << bodyName << endl;
    }
    const dictionary& bodyStateDict = floaterStateDict.subDict(bodyName);

    // Reading motionSolver coeff dict from dynamicMeshDict needed for creating
    // floaterMotion object.
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    const dictionary& solverDict =
        dynamicMeshDict.subDict("solvers").subDict(bodyName);

    word motionSolverType = solverDict.get<word>("motionSolver");
    const dictionary& motionSolverCoeffDict =
        solverDict.subDict(motionSolverType + "Coeffs");

    // Creating floaterMotion object to hold rigid body state and calculate
    // accelerations.
    floaterMotion body(motionSolverCoeffDict, bodyStateDict, runTime);

    // Reading active degrees of freedom from motionSolverCoeffDict (copied from
    // floaterMeshMotionSolver constructor)
    Vector<bool> linDirs(motionSolverCoeffDict.getOrDefault
    (
        "linDirs", Vector<bool>(1, 1, 1))
    );
    Vector<bool> rotDirs(motionSolverCoeffDict.getOrDefault
    (
        "rotDirs", Vector<bool>(1, 1, 1))
    );

    // Ensure no linear motion along and no rotational motion transverse to
    // empty directions
    const Vector<label> solDirs = (mesh.solutionD() + Vector<label>::one)/2;
    forAll(solDirs, i)
    {
        if (!solDirs[i])
        {
            linDirs[i] = 0;
            rotDirs[(i + 1) % 3] = 0;
            rotDirs[(i + 2) % 3] = 0;
        }
    }

    // Adding linear degrees of freedom to DoFs
    labelList DoFs;
    forAll(linDirs, i)
    {
        if (linDirs[i])
        {
            DoFs.append(i);
        }
    }

    // Adding rotational degrees of freedom to DoFs
    forAll(rotDirs, i)
    {
        if (rotDirs[i])
        {
            DoFs.append(i+3);
        }
    }
   Info << "Active degrees of freedom: " << DoFs << endl;

    // Creating ODE system for rigid body.
    kelvinKirchhoffEquations ode(body, DoFs);

    // Create the selected ODE system solver
    dictionary ODESolverDict;
    ODESolverDict.add("solver", "RKDP45"); //Note: hardcoded ODEsolver
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, ODESolverDict);

    odeSolver->relTol() = 1e-12; //Note: Hardcoded tolerances
    odeSolver->absTol() = 1e-12;
    scalar dtEst = 1e-8;

    // Initialise the ODE system fields
    scalarField Y(ode.nEqns(),0);
    // State: Y = [x0, Q, v0, omega]
    // nDim = 3+9+3+3=18

    for (label n = 0; n<3; n++)
    {
        Y[n] = body.centreOfRotation()[n];
        Y[n+12] = body.v()[n];
        Y[n+12+3] = body.omega()[n];
    }
    label k = 3;
    tensor Q = body.orientation();
    for (label m = 0; m<3; m++)
    {
        for (label n = 0; n<3; n++)
        {
            Y[k] = Q(m,n);
            k++;
        }
    }

    // Making file in postProcessing to write body state at every time step
    fileName filePath = runTime.path() / "postProcessing" / "floaterMotion"
        / runTime.timeName() / "floaterStateWriter.dat";
    mkDir(filePath.path());
    OFstream os(filePath);
    os << "# Rigid body motion state for body" << endl;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        odeSolver->solve(0, runTime.deltaTValue(), Y, dtEst);
        Info << "Solution: t = " << runTime.timeName() << ", Y = " << Y << endl;

        // Converting from Y to x0, Q, v0, omega and inserting in body object
        vector x0(Y[0], Y[1], Y[2]);
        tensor Q
        (
            Y[3], Y[4], Y[5],
            Y[6], Y[7], Y[8],
            Y[9], Y[10], Y[11]
        );
        vector v0(Y[12], Y[13], Y[14]);
        vector omega(Y[15], Y[16], Y[17]);
        body.setCentreOfRotation(x0);
        body.setOrientation(Q);
        body.setVelocity(v0);
        body.setAngularVelocity(omega);

        // Writing body state to postProcessing file
        const floaterMotionState& state = body.state();
        os  << "Time = " << runTime.timeName() << endl;
        os  << "Centre of rotation  : " << state.centreOfRotation() << endl
            << "Velocity            : " << state.v() << endl
            << "Acceleration        : " << state.a() << endl
            << "Orientation         : " << state.Q() << endl
            << "Angular velocity    : " << state.omega() << endl
            << "Angular acceleration: " << state.domegadt() << endl
            << "Centre of mass      : " << body.centreOfMass() << endl
            << "CoM velocity        : " << body.velocity(body.centreOfMass()) << endl
            << endl;

        // Write floaters file to time directory
        if (runTime.writeTime())
        {
            IOdictionary dict
            (
                IOobject
                (
                    "floaters",
                    mesh.time().timeName(),
                    "uniform",
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

            // Create or get the sub-dictionary with the name from motion_.name()
            dictionary& bodyOutDict = dict.subDictOrAdd(body.name());
            body.write(bodyOutDict);

            dict.regIOobject::write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
