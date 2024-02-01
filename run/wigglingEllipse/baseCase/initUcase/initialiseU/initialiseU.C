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
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

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

//    while (simple.loop())
//    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Solve Laplace equation for velocity potential 
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn(fvm::laplacian(Phi));

            PhiEqn.solve();
        }

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
//    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
