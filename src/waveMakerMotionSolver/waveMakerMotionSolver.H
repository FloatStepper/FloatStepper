/*---------------------------------------------------------------------------*\
    Module Name:       FloatStepper
    Description:       OpenFOAM extension module for fluid-rigid body coupling
    License:           GNU General Public License (GPL) version 3
    Copyright:         2023 Johan Roenby, STROMNING APS
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
    Copyright (C) 2023 Johan Roenby
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

Class
    Foam::waveMakerMotionSolver

Group
   grpMeshMotionSolvers

Description
    Wave maker mesh motion solver for an fvMesh.

SourceFiles
    waveMakerMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef waveMakerMotionSolver_H
#define waveMakerMotionSolver_H

#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class waveMakerMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class waveMakerMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Array of times
        scalarField times_;

        //- Array of piston positions corresponding to times_
        List<scalarField> pistonPositions_;
        
        //- Number of pistons
        int nPistons_;

        //- Left bound of region with extracting/contracting cells
        scalar xl_;
        
        //- Right bound of region with extracting/contracting cells
        scalar xr_;

        //- Number of times to repeat tabulated piston motion (default = 1) 
        label repetitions_;

        //- Interpolation method for piston positions between listed times
        //  (default = spline) 
        word timeInterpolation_;

        //- Interpolation method for positions between listed piston positions
        //  (default = linear) 
        word spaceInterpolation_;

        //- Y coordinate of piston centres
        scalarField yPistonCentres_;

        //- List of z positions for scaling motion e.g. to generate hinged wave
        //  maker motion.
        scalarField zPositions_;

        //- Degree of scaling of motion at z positions (linear interpolation in
        //  between).
        scalarField zScaling_;

        //- Number of moving mesh points 
        int nMovingPoints_;    
        
        //- Indices of moving points in points file
        labelList movingPoints_;
        
        //- Original x-coordinates of moving points
        scalarField xOriginal_;
        
        //- y-coordinate of moving points
        scalarField yOriginal_;

        //- Scaling as function of height (calculated from zPositions and
        //  zScaling).
        scalarField scaling_;

        //- Should piston positions be written to log (default: no) 
        bool writePositionsToLogFile_;

        //- Optional factor multiplied on piston positions to amplify or
        //  attenuate them
        scalar amplificationFactor_;


    // Private Member Functions

        //- No copy construct
        waveMakerMotionSolver
        (
            const waveMakerMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const waveMakerMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("waveMakerMotionSolver");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        waveMakerMotionSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~waveMakerMotionSolver() = default;


    // Member Functions

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
