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
    Copyright (C) 2019 OpenCFD Ltd.
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
    floaterMotionPotentialFlow

Group
    grpBasicSolvers

Description
    Laplace equation solver.

    \heading Solver details
    The solver is applicable to, e.g. the velocity potential of a moving rigid
    body.  The equation is given by:

    \f[
        \div \left( \grad phi \right) = 0
    \f]

    Where:
    \vartable
        phi     | is the velocity potential.
    \endvartable

    \heading Required fields
    \plaintable
        phi     | Scalar field which is solved for, e.g. velocity potential
    \endplaintable

    Based on laplacianFoam.

    Johan Roenby, 2025

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
//#include "rigidBodyVelocityPotentialFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar velocity potential."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating potential flow for rigid body motion...\n" << endl;

    // Reading dynamicMeshDict
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
    
    // Iterate through all bodies
    forAll(bodies, i)
    {
        const word& body = bodies[i];

        // Ignore entries that are not subdicts
        if (!floaterStateDict.isDict(body)) continue;

        // Access the sub-dictionary
        const dictionary& bodyDict = floaterStateDict.subDict(body);

        // Ignore subDicts not containing relevant body info
        if
        (
            !bodyDict.found("centreOfRotation") 
         || !bodyDict.found("velocity")
         || !bodyDict.found("angularVelocity")
        ) continue;

        Info << "Reading body velocity for body names: " << body << endl;

        // Reading body position and velocity from floaters file
        point CoR = bodyDict.get<point>("centreOfRotation");
        vector v0 = bodyDict.get<vector>("velocity");
        vector omega = bodyDict.get<vector>("angularVelocity");
        
        // Reading body patches
        const dictionary& bodySubDict =
            dynamicMeshDict.subDict("solvers").subDict(body);
        const word motionSolverType = bodySubDict.get<word>("motionSolver");
        const dictionary& solverDict =
            bodySubDict.subDict(motionSolverType + "Coeffs");
        wordRes patches = solverDict.get<wordRes>("patches");

        // Setting velocity potential boundary condition on floater
        const auto& bMesh = mesh.boundary();
        forAll(patches, patchi)
        {
            const label patchId
                = mesh.boundaryMesh().findPatchID(patches[patchi], false);

            // Access the boundary condition and cast to fixedGradient
            fixedGradientFvPatchField<scalar>& Phipatch =
            refCast<fixedGradientFvPatchField<scalar>>(Phi.boundaryFieldRef()[patchId]);
            // Set the gradient (e.g., to 100)
            const vectorField nf(bMesh[patchId].nf());
            const vectorField Cf(bMesh[patchId].Cf());
            Phipatch.gradient() = (nf & (v0 + (omega ^ (Cf - CoR))));
        }
    }

    // Solve Laplace equation for velocity potential 
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix PhiEqn(fvm::laplacian(Phi));

        PhiEqn.solve();
    }

    // Calculate velocity field from potential
    volVectorField gradPhi(fvc::grad(Phi));
    const vectorField& gradPhiIn = gradPhi.primitiveField();
    vectorField& UIn = U.primitiveFieldRef();
    forAll(UIn, celli)
    {
        UIn[celli] = gradPhiIn[celli];
    }
//        U.correctBoundaryConditions();
    U.write();

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
