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
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
    Copyright (C) 2023 Johan Roenby
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

#include "waveMakerMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateSplineXY.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(waveMakerMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        waveMakerMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveMakerMotionSolver::waveMakerMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    times_(coeffDict().lookup("times")),
    pistonPositions_(coeffDict().lookup("pistonPositions")),
    nPistons_(pistonPositions_.size()),
    xl_(readScalar(coeffDict().lookup("xLeft"))),
    xr_(readScalar(coeffDict().lookup("xRight"))),
    repetitions_(coeffDict().lookupOrDefault("repetitions", 0)),
    timeInterpolation_
    (
        coeffDict().lookupOrDefault<word>("timeInterpolation", "spline")
    ),
    spaceInterpolation_
    (
        coeffDict().lookupOrDefault<word>("spaceInterpolation", "linear")
    ),
    yPistonCentres_(nPistons_),
    zPositions_(coeffDict().lookup("zPositions")),
    zScaling_(coeffDict().lookup("zScaling")),
    nMovingPoints_(sum(neg(points0().component(vector::X) - xr_))),
    movingPoints_(nMovingPoints_),
    xOriginal_(nMovingPoints_),
    yOriginal_(nMovingPoints_),
    scaling_(nMovingPoints_),
    writePositionsToLogFile_
    (
        coeffDict().lookupOrDefault<bool>("positionsToLog", 0)
    ),
    amplificationFactor_(coeffDict().lookupOrDefault("amplification", 1.0))
{

    //Amplifying/attenuating signal
    pistonPositions_ = amplificationFactor_*pistonPositions_;
    
    //Finding moving mesh points
    const pointField& p = points0();
    scalarField isMoving(neg(p.component(vector::X) - xr_));
    label iMov(-1);
    forAll(isMoving, ip)
    {
        if (isMoving[ip])
        {
            iMov++;
            movingPoints_[iMov] = ip;
            xOriginal_[iMov] = p[ip].component(vector::X);
            yOriginal_[iMov] = p[ip].component(vector::Y);
        }
    }

    scaling_ = 1;
    if (zPositions_.size() > 1)
    {
        forAll(movingPoints_, ip)
        {
            scalar zp = p[movingPoints_[ip]].component(vector::Z);
            scaling_[ip] = Foam::interpolateXY(zp, zPositions_, zScaling_);
        }
    }

    //Reading piston y-positions
    if (coeffDict().found("yPiston")) 
    {
        yPistonCentres_ = readScalar(coeffDict().lookup("yPiston"));
    } 
    else
    {
        Info << "yPiston not found. Distributing pistons uniformly between y = " 
            << gMax(yOriginal_) << " and "
            << gMin(yOriginal_) << endl;

        scalar pistonWidth = (gMax(yOriginal_)-gMin(yOriginal_))/nPistons_ ;
        forAll(yPistonCentres_,ip)
        {
            yPistonCentres_[ip] = gMin(yOriginal_) + (ip+0.5)*pistonWidth;
        }
    }
    if (writePositionsToLogFile_)
    {
        Info << "Piston y-positions: " << yPistonCentres_ << endl;
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::pointField>
Foam::waveMakerMotionSolver::curPoints() const
{
    tmp<pointField> newPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    if (!moveAllCells())
    {
        Info << "Note: not moving all cells - returning transformed points."
            << endl;
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs()) =
            pointField(newPoints.ref(), pointIDs());

        return ttransformedPts;
    }

    return newPoints;
}

void Foam::waveMakerMotionSolver::solve()
{

    // For times later than latest time in tabulated table the mesh is treated
    // as static unless a number of repetitions is specified
    scalar t = time().value();
    const scalar tmin(min(times_));
    const scalar tmax(max(1, repetitions_)*max(times_));
    t = max(tmin, min(t, tmax));

    if (repetitions_ > 0)
    {
        const scalar period = tmax/repetitions_;
        t = fmod(t, period);
    }

    // Reference to point displacements
    vectorField& pd = pointDisplacement_.primitiveFieldRef();

    Info << "Moving points in waveMakeMotionSolver" << endl;
    //Finding piston positions xPist at current time by interpolation
    scalarField xPist(nPistons_);
    forAll(xPist, ip)
    {
        if ( timeInterpolation_ == "spline" )
        {
            xPist[ip] =
                Foam::interpolateSplineXY(t, times_, pistonPositions_[ip]);
        }
        else if (timeInterpolation_ == "linear" )
        {
            xPist[ip] =
                Foam::interpolateXY(t, times_, pistonPositions_[ip]);
        }
        else
        {
            Info << "Warning: unknown timeInterpolation method (options: "
                << "linear and spline)" << endl;                
        }
    }
    
    if (writePositionsToLogFile_)
    {
        Info << "Piston x-positions: " << xPist << endl;
    }

    // Creating piston position x and y arrays for use in subsequent
    // interpolation.
    scalarField Xi(nPistons_+2), Yi(nPistons_+2);
    Xi[0] = xPist[0];
    Yi[0] = -1e10;
    forAll(xPist, ix) 
    {
        Xi[ix+1] = xPist[ix];
        Yi[ix+1] = yPistonCentres_[ix];
    }
    Xi[nPistons_+1] = xPist[nPistons_-1];
    Yi[nPistons_+1] = 1e10;

    scalarField XNEW(nMovingPoints_);
    if ( spaceInterpolation_ == "spline" )
    {
        XNEW = Foam::interpolateSplineXY(yOriginal_, Yi, Xi);
    }
    else if ( spaceInterpolation_ == "linear" )
    {
        XNEW = Foam::interpolateXY(yOriginal_, Yi, Xi);
    }
    else
    {
        Info << "Warning: unknown spaceInterpolation method (options: "
            << "linear and spline)" << endl;
    }

    //Changing x-coordinate of all moving points
    scalar xmid(0.5*(xl_+xr_)), L(0.5*(xr_-xl_)), X(0.0);
    forAll(movingPoints_, ip)
    {    
        scalar x = xOriginal_[ip];
        scalar xOrig = x;
        if (x < xl_) 
        {
            x = (x + XNEW[ip]);
        }
        else if ( x >= xl_ && x < xmid ) 
        {
            X = x - xl_;
            X = X*(1 - .5*(X/L)*(XNEW[ip]/L));
            x = X + xl_ + XNEW[ip];
        }
        else if ( x >= xmid && x < xr_ ) 
        {
            X = -(x - xr_);
            X = X*(1 - 0.5*(X/L)*(XNEW[ip]/L));
            x = xr_ - X;
        }
        scalar dx = scaling_[ip]*(x - xOrig);
        pd[movingPoints_[ip]].replace(vector::X, dx);
    }
}


bool Foam::waveMakerMotionSolver::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
/*
    IOdictionary dict
    (
        IOobject
        (
            "waveMakerMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    motion_.state().write(dict);
    return dict.regIOobject::write();
*/
    return false;
}


bool Foam::waveMakerMotionSolver::read()
{
/*
    if (displacementMotionSolver::read())
    {
        motion_.read(coeffDict());

        return true;
    }
*/
    return false;
}


// ************************************************************************* //
