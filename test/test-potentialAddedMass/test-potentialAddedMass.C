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

Application
    test-potentialAddedMass

Group

Description
    Testing added mass calculator based on velocity potential.

Author
    Johan Roenby & Henning Scheufler

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "addedMass.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "calculates the added mass"
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating added mass matrix\n" << endl;

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );


    autoPtr<addedMass> addedM = addedMass::New
    (
        dynamicMeshDict.get<word>("addedMassModel"),
        mesh,
        dynamicMeshDict
    );

    // Determining active degrees of freedom
    
    Vector<bool> linDirs(dynamicMeshDict.getOrDefault
    (
        "linDirs", Vector<bool>(1, 1, 1))
    );
    Vector<bool> rotDirs(dynamicMeshDict.getOrDefault
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

    // Reading center of rotation and orientation matrix from dicitonary
    vector CoR = dynamicMeshDict.get<vector>("CofR");
    tensor Q = dynamicMeshDict.get<tensor>("Q");
    wordRes patches = dynamicMeshDict.get<wordRes>("patches");

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
