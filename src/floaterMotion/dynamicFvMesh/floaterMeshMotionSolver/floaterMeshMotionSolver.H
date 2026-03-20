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
    Foam::floaterMeshMotionSolver

Group
   grpMeshMotionSolvers

Description
    Description

    Motion solver for floating body using six degrees of freedom rigid
    body motion with customised box-shaped morphing regions for the linear
    displacement and SLERP for the rotational morphing.

    Between the inner and outer boxes the mesh morphing transitions linearly
    from full rigid body motion to stationary. Inside the inner box the mesh
    moves as a rigid body. Outside the outer box the mesh is stationary.

    Activated by setting separateMeshMorphing true. If false, falls back to
    standard distance-based blending using innerDistance and outerDistance.
    The y and z directions inherit all settings from xMorphing if not
    specified.

    Usage
    \verbatim
    solvers
    {
        box
        {
            motionSolverLibs    ( floaterMotion );
            motionSolver        floaterMeshMotionSolver;

            floaterMeshMotionSolverCoeffs
            {
                patches                 ( box );
                innerDistance           0.075;
                outerDistance           0.4;
                separateMeshMorphing    true;

                xMorphing
                {
                    coordinateSystem { origin (0 0 0); e1 (1 0 0); e2 (0 1 0); }
                    relativeTo       body;
                    innerBox     { xMin -0.6; xMax 0.5; zMin -0.21; zMax 0.19; }
                    outerBox     { xMin -2.2; xMax 2.2; zMin -0.39; zMax 0.39; }
                }

                yMorphing { $xMorphing; }

                zMorphing
                {
                    $xMorphing;
                    innerBox     { xMin -0.6; xMax 0.5; zMin -0.21; zMax 0.19; }
                    outerBox     { xMin -2.2; xMax 2.2; zMin -0.39; zMax 0.39; }
                }
            }
        }
    }
    \endverbatim

    Note
    Bounds not specified in innerBox and outerBox default to +/-GREAT and
    +/-VGREAT, respectively, ensuring outer bounds always exceed inner bounds
    for unspecified directions. The coordinateSystem entry defaults to the
    global Cartesian system if omitted. relativeTo body means bounds move
    with the body; relativeTo lab means bounds are fixed in the global frame.

SourceFiles
    floaterMeshMotionSolver.C, morphingBoxRegion.C

\*---------------------------------------------------------------------------*/

#ifndef floaterMeshMotionSolver_H
#define floaterMeshMotionSolver_H

#include "pointDisplacementMotionSolver.H"
#include "floaterMotion.H"
#include "morphingBoxRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class floaterMeshMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class floaterMeshMotionSolver
:
    public pointDisplacementMotionSolver
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

        //- Switch to separate linear and rotational mesh deformation treatment
        const bool separateMeshMorphing_;

        //- First horisontal direction in mesh motion coordinates used if
        //  translation and rotation mesh deformation are treated separately
        vector xhat_;

        //- Second horisontal direction in mesh motion coordinates used if
        //  translation and rotation mesh deformation are treated separately
        vector yhat_;

        //- Vertical direction in mesh motion coordinates used if translation
        //  and rotation mesh deformation are treated separately
        vector zhat_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        //scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Mesh morphing blending between stationary and co-moving
        pointScalarField scale_;

        //- Mesh morphing blending between stationary and co-moving for linear
        //  displacement. Only allocated if separateMeshMorphing_.
        autoPtr<pointScalarField> xScale_;
        autoPtr<pointScalarField> yScale_;
        autoPtr<pointScalarField> zScale_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- List of active DoFs (e.g. [0 1 5] for 2D motion with rotation)
        labelList DoFs_;

    // Private Member Functions

        //- No copy construct
        floaterMeshMotionSolver
        (
            const floaterMeshMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const floaterMeshMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("floaterMeshMotionSolver");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        floaterMeshMotionSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~floaterMeshMotionSolver() = default;


    // Member Functions

        //- Return the six DoF motion object
        const floaterMotion& motion() const;

        //- Return the six DoF motion object
        floaterMotion& motion();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update pointDisplacement field based on floater position and
        //  orientation
        void updateDisplacement();

        //- Apply grading to stretched morphing zone
        void applyGrading
        (
            const vector& dir,
            const scalarField& pts,
            const scalar Lx,
            scalarField& scale
        ) const;

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
