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
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "dynamicPointDisplacementFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "pointMesh.H"
#include "volFields.H"
#include "pointDisplacementMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicPointDisplacementFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicPointDisplacementFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicPointDisplacementFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicPointDisplacementFvMesh::dynamicPointDisplacementFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit),
    points0_
    (
        IOobject
        (
            "points0",
            time().constant(),
            polyMesh::meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointIOField
        (
            IOobject
            (
                "points",
                time().constant(),
                polyMesh::meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    ),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement",
            time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(*this)
    ),
    motionSolvers_(),
    onlyMeshMotion_(false)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicPointDisplacementFvMesh::init
(
    const bool doInit,
    const bool mandatory
)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    IOobject ioDict
    (
        "dynamicMeshDict",
        time().constant(),
        *this,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    IOdictionary dict(ioDict);

    // Read onlyMeshMotion used to only move mesh with moveDynamicMesh
    onlyMeshMotion_ = dict.getOrDefault<bool>("onlyMeshMotion", false);

    label i = 0;
    if (dict.found("solvers"))
    {
        const dictionary& solvertDict = dict.subDict("solvers");


        motionSolvers_.setSize(solvertDict.size());

        for (const entry& dEntry : solvertDict)
        {
            if (dEntry.isDict())
            {
                IOobject io(ioDict);
                io.readOpt(IOobject::NO_READ);
                io.writeOpt(IOobject::AUTO_WRITE);
                io.rename(dEntry.dict().dictName());

                IOdictionary IOsolverDict
                (
                    io,
                    dEntry.dict()
                );

                motionSolvers_.set
                (
                    i++,
                    motionSolver::New(*this, IOsolverDict)
                );
            }
        }
        motionSolvers_.setSize(i);
    }
    else if (mandatory)
    {
        motionSolvers_.setSize(1);
        motionSolvers_.set(i++, motionSolver::New(*this, dict));
    }

    // Assume something changed
    return true;
}


bool Foam::dynamicPointDisplacementFvMesh::init(const bool doInit)
{
    return init(doInit, true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicPointDisplacementFvMesh::~dynamicPointDisplacementFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dynamicPointDisplacementFvMesh::mapFields
(
    const mapPolyMesh& mpm
)
{
    dynamicFvMesh::mapFields(mpm);

    // Update the motionSolvers for any topo change ...
    for (auto& ms : motionSolvers_)
    {
        ms.updateMesh(mpm);
    }
}


bool Foam::dynamicPointDisplacementFvMesh::update()
{
    if (motionSolvers_.size())
    {
        // Accumulated displacement
        pointField& disp = pointDisplacement_;
        disp = Zero;

        // Accumulate all point displacements
        forAll(motionSolvers_, i)
        {
            if (isA<pointDisplacementMotionSolver>(motionSolvers_[i]))
            {
                pointDisplacementMotionSolver& pointDistSolver =
                    dynamic_cast<pointDisplacementMotionSolver&>
                    (
                        motionSolvers_[i]
                    );
                pointDistSolver.solve();
            }
        }

        // Moving mesh with accummulated point displacement, disp
        fvMesh::movePoints(points0_ + disp);

        // Note: pointDisplacement_ boundary conditions are needed for restart
        // but internal values are recalculated. Therefore setting to zero
        // before write out to save disk space.
        disp = Zero;

        volVectorField* Uptr = getObjectPtr<volVectorField>("U");

        if (Uptr)
        {
            Uptr->correctBoundaryConditions();
        }
    }

    return true;
}


// ************************************************************************* //