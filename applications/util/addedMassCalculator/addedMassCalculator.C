/*--------------------------------------------------- -----------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*--------------------------------------------------- -----------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
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
    addedMassCalculator


Description
    Application for calculating the added mass of a rigid body.
    Based on interFoam to get its createFields.H.

Author
    Johan Roenby, 2025

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "addedMass.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculates the added mass of a body whose patches are specified in\n"
        "dynamicMeshDict under solvers.<bodyName>.\n"
        "<bodyName> must be specified using the -body argument."
    );

    argList::addOption
    (
        "body",
        "body",
        "Specify the name of the body whose added mass to calculate."
        " (no default)"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"

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
    const word body = args.get<word>("body");
    const dictionary& bodyDict =
        dynamicMeshDict.subDict("solvers").subDict(body);
    const word motionSolverType = bodyDict.get<word>("motionSolver");
    const dictionary& solverDict =
        bodyDict.subDict(motionSolverType + "Coeffs");

    autoPtr<addedMass> addedM = addedMass::New
    (
        solverDict.get<word>("addedMassModel"),
        mesh,
        dynamicMeshDict
    );

    wordRes patches = solverDict.get<wordRes>("patches");

    // Determining active degrees of freedom
    
    Vector<bool> linDirs(solverDict.getOrDefault
    (
        "linDirs", Vector<bool>(1, 1, 1))
    );
    Vector<bool> rotDirs(solverDict.getOrDefault
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

    FixedList<bool,6> enabledDoFs(false);

    // Adding linear degrees of freedom to DoFs_
    forAll(linDirs, i)
    {
        if (linDirs[i])
        {
            enabledDoFs[i] = true;
        }
    }

    // Adding rotational degrees of freedom to DoFs_
    forAll(rotDirs, i)
    {
        if (rotDirs[i])
        {
            enabledDoFs[i+3] = true;
        }
    }
   Info << "Active degrees of freedom: " << enabledDoFs << endl;

    OFstream osI("ma.dat");

    // Read floaterMotionState
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

    vector CoR = floaterStateDict.subDict(body).get<vector>("centreOfRotation");
    tensor Q = floaterStateDict.subDict(body).get<tensor>("orientation");

    // Added mass matrix
    typedef SquareMatrix<scalar> SMatrix;
    SMatrix Madd({6, 6}, 0);

    Madd = addedM->computeAddedMass
    (
        enabledDoFs,
        CoR,
        Q,
        patches
    );

    osI << " " << Madd << endl;
    Info << "Madd = " << Madd << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
