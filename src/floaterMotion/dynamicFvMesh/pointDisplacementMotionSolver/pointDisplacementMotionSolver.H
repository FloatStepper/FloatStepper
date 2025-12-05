/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

Class
    Foam::pointDisplacementMotionSolver

Description
    Virtual base class for displacement motion solvers, where the point
    motion is relative to a set of fixed points (points0).

    Differs from points0MotionSolver by all instances sharing a reference to the
    same points0 pointIOField and same pointDisplacement pointVectorField.
    These fields are read from object registry on construction so either solver
    or dynamicFvMesh class must have created them.

SourceFiles
    pointDisplacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef pointDisplacementMotionSolver_H
#define pointDisplacementMotionSolver_H

#include "motionSolver.H"
#include "pointFields.H"

namespace Foam
{

class pointDisplacementMotionSolver
:
    public motionSolver
{
protected:

    //- Reference to the points0 field (owned by registry)
    const pointIOField& points0_;

    //- The specific displacement field calculated by this solver instance
    pointVectorField& pointDisplacement_;

private:

    //- Helper: Check registry for points0, or load from constant/polyMesh/points
    static const pointIOField& lookupOrLoadPoints0(const polyMesh& mesh);

    //- Helper: Generate unique name based on dictionary name (e.g., dispField_box1)
    static word determineFieldName(const IOdictionary& dict);

public:

    //- Runtime type information
    TypeName("pointDisplacementMotionSolver");

    // Constructors

    //- Construct from mesh and dictionary
    pointDisplacementMotionSolver
    (
        const polyMesh&,
        const IOdictionary&,
        const word& type
    );

    //- Destructor
    virtual ~pointDisplacementMotionSolver();

    // Member Functions

    //- Return reference to the internal displacement field
    virtual pointVectorField& pointDisplacement()
    {
        return pointDisplacement_;
    }

    //- Return point location: points0 + pointDisplacement_
    virtual tmp<pointField> curPoints() const;

    //- Solve for motion
    virtual void solve();

    //- Update topology
    virtual void updateMesh(const mapPolyMesh&);

    //- Update map (Mandatory override from base class)
    virtual void movePoints(const pointField&);
};

} // End namespace Foam

#endif
