/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "pointDisplacementMotionSolver.H"
//#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointDisplacementMotionSolver, 0);
//    addToRunTimeSelectionTable(motionSolver, pointDisplacementMotionSolver, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pointIOField& Foam::pointDisplacementMotionSolver::lookupOrLoadPoints0
(
    const polyMesh& mesh
)
{
    // 1. Check if points0 already exists in the objectRegistry
    if (mesh.foundObject<pointIOField>("points0"))
    {
        DebugInfo << "pointDisplacementMotionSolver: Reusing existing points0 from registry." << endl;
        return mesh.lookupObject<pointIOField>("points0");
    }

    // 2. Load from constant/polyMesh/points explicitly
    // We strictly use mesh.time().constant() to avoid reading from 0/ or other time dirs.
    IOobject io
    (
        "points",                   // Look for file named "points"
        mesh.time().constant(),     // STRICTLY in constant directory
        polyMesh::meshSubDir,       // inside polyMesh/
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE          // Never overwrite the constant points
    );

    pointIOField* p0Ptr = nullptr;

    if (io.typeHeaderOk<pointIOField>(true))
    {
        Info<< "pointDisplacementMotionSolver: Loading points0 from " 
            << io.path() << endl;
            
        p0Ptr = new pointIOField(io);
        
        // Rename the object in memory to "points0" so it doesn't clash with "points"
        p0Ptr->rename("points0");
    }
    else
    {
        // Fallback: If constant/polyMesh/points doesn't exist (rare), 
        // use current mesh points.
        Info<< "pointDisplacementMotionSolver: constant/polyMesh/points not found. "
            << "Using current mesh points as points0." << endl;
            
        p0Ptr = new pointIOField
        (
            IOobject
            (
                "points0",
                mesh.time().constant(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh.points()
        );
    }

    // Transfer ownership to the registry so other solvers can access it
    p0Ptr->store();

    return *p0Ptr;
}

Foam::word Foam::pointDisplacementMotionSolver::determineFieldName
(
    const IOdictionary& dict
)
{
    // Try to get the dictionary name (e.g., "box1" from the solvers dict)
    word dictName = dict.dictName();

    // If dictName is empty or generic, look for an explicit "name" keyword
    // If that fails, default to "pointDisplacement"
    if (dictName.empty() || dictName == "solvers")
    {
        return dict.lookupOrDefault<word>("name", "pointDisplacement");
    }
    
    // Automatic naming: dispField_box1
    return "dispField_" + dictName;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointDisplacementMotionSolver::pointDisplacementMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    points0_(lookupOrLoadPoints0(mesh)),
    pointDisplacement_
    (
        mesh.lookupObjectRef<pointVectorField>("pointDisplacement")
    )
{
    // Initialize if not read from disk (e.g. starting from 0)
    if (!pointDisplacement_.headerOk())
    {
        pointDisplacement_ == vector::zero;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointDisplacementMotionSolver::~pointDisplacementMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::pointField> Foam::pointDisplacementMotionSolver::curPoints() const
{
    // Return Original Points + Local Displacement
    tmp<pointField> tnewPoints(new pointField(points0_));
    tnewPoints.ref() += pointDisplacement_.internalField();
    
    return tnewPoints;
}


void Foam::pointDisplacementMotionSolver::solve()
{}


void Foam::pointDisplacementMotionSolver::movePoints(const pointField&)
{}


void Foam::pointDisplacementMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    // 1. Update the mesh and automatically map the registered pointDisplacement_ field
    motionSolver::updateMesh(mpm);

    // 2. Handle points0
    // Since points0 is shared, multiple solvers will call this function.
    // We must ensure points0 is only updated ONCE.
    // We do this by checking if the size of points0 matches the NEW mesh size.
    
    // We cast away const because we are responsible for maintaining this shared object
    pointIOField& p0 = const_cast<pointIOField&>(points0_);
    
    // If sizes match, a previous solver (or this one) has already updated it.
    if (p0.size() == mpm.pointMap().size())
    {
        return;
    }

    Info<< "pointDisplacementMotionSolver: Updating shared points0 field for topology change." << endl;

    // Adapted from points0MotionSolver

    // Get the new points either from the map or the mesh
    const pointField& points =
    (
        mpm.hasMotionPoints()
      ? mpm.preMotionPoints()
      : mesh().points()
    );

    // Note: boundBox does reduce
    const vector span0 = boundBox(p0).span();
    const vector span = boundBox(points).span();

    vector scaleFactors(cmptDivide(span0, span));

    pointField newPoints0(mpm.pointMap().size());

    forAll(newPoints0, pointi)
    {
        label oldPointi = mpm.pointMap()[pointi];

        if (oldPointi >= 0)
        {
            label masterPointi = mpm.reversePointMap()[oldPointi];

            if (masterPointi == pointi)
            {
                newPoints0[pointi] = p0[oldPointi];
            }
            else
            {
                // New point - assume motion is scaling
                // This logic places the new point in the points0 configuration
                // based on its relative position in the current mesh
                newPoints0[pointi] = p0[oldPointi] + cmptMultiply
                (
                    scaleFactors,
                    points[pointi] - points[masterPointi]
                );
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot determine coordinates of introduced vertices."
                << " New vertex " << pointi << " at coordinate "
                << points[pointi] << exit(FatalError);
        }
    }

    // Apply 2D corrections if strictly 2D case
    twoDCorrectPoints(newPoints0);

    // Transfer data to the shared object
    p0.transfer(newPoints0);

    // Update metadata for the shared object
    p0.rename("points0");
    p0.writeOpt(IOobject::AUTO_WRITE);
    // It is important to set the instance to current time so it follows the mesh
    // into the new time directory
    p0.instance() = time().timeName(); 
}

/*
void Foam::pointDisplacementMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    motionSolver::updateMesh(mpm);
}
*/

// ************************************************************************* //
