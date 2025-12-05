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
    Copyright (C) 2020,2022 OpenCFD Ltd.
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

Class
    Foam::dynamicPointDisplacementFvMesh

Description
    Dynamic mesh able to handle multiple motion solvers.
    NOTE: If the word entry "solvers" is not found it falls back to a single
    motion solver behavior.

SourceFiles
    dynamicPointDisplacementFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicPointDisplacementFvMesh_H
#define dynamicPointDisplacementFvMesh_H

#include "dynamicFvMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                Class dynamicPointDisplacementFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicPointDisplacementFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- The base points from which point displacements are calculated
        const pointIOField points0_;

        //- Field to hold accumulated mesh point displacement
        pointVectorField pointDisplacement_;

        // List of motionSolvers
        PtrList<motionSolver> motionSolvers_;

        // Switch to allow initial pure mesh motion with moveDynamicMesh
        bool onlyMeshMotion_;

    // Private Member Functions

        //- No copy construct
        dynamicPointDisplacementFvMesh
        (
            const dynamicPointDisplacementFvMesh&
        ) = delete;

        //- No copy assignment
        void operator=(const dynamicPointDisplacementFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicPointDisplacementFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicPointDisplacementFvMesh
        (
            const IOobject& io,
            const bool doInit=true
        );


    //- Destructor
    ~dynamicPointDisplacementFvMesh();


    // Member Functions

        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        //- Initialise all non-demand-driven data. Used for initialising
        //  without default motionSolver (used by dynamicRefineFvMesh)
        virtual bool init(const bool doInit, const bool mandatory);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Map all fields in time using given map. Triggered by topo change
        virtual void mapFields(const mapPolyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //