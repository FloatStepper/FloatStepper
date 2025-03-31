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
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::floaterMotionSolver

Group
   grpMeshMotionSolvers

Description
    6-DoF solid-body mesh motion solver for an fvMesh.

    Applies SLERP interpolation of movement as function of distance to
    the object surface.

SourceFiles
    floaterMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef floaterMotionSolver_H
#define floaterMotionSolver_H

#include "displacementMotionSolver.H"
#include "floaterMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class floaterMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class floaterMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Six DoF motion object
        floaterMotion motion_;

        //- Names of patches making up the rigid body
        wordRes patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Optional stretching of distance-to-floater measure to accomodate
        //  different effective di_ and do_ in different directions
        const tensor distStretch_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        //scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- List of active DoFs (e.g. [0 1 5] for 2D motion with rotation)
        labelList DoFs_;

    // Private Member Functions

        //- No copy construct
        floaterMotionSolver
        (
            const floaterMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const floaterMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("floaterMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        floaterMotionSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~floaterMotionSolver() = default;


    // Member Functions

        //- Return the six DoF motion object
        const floaterMotion& motion() const;

        //- Return the six DoF motion object
        floaterMotion& motion();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();


        //- Write state using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;


        //- Read dynamicMeshDict dictionary
        virtual bool read();


        //- Return list of active DoFs
        const labelList& DoFs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
