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

\*---------------------------------------------------------------------------*/

#include "floaterMeshMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(floaterMeshMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        floaterMeshMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMeshMotionSolver::floaterMeshMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    pointDisplacementMotionSolver(mesh, dict, typeName),
    motion_
    (
        coeffDict(), //Inherited from motionSolver - reads from dynamicMeshDict
        IOdictionary
        (
            IOobject
            (
                "floaters",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(dict.dictName()),
        mesh.time()
    ),
    patches_(coeffDict().get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(coeffDict().get<scalar>("innerDistance")),
    do_(coeffDict().get<scalar>("outerDistance")),
    distStretch_(coeffDict().getOrDefault<tensor>("distStretch", tensor::I)),
    separateTransRotPointDisplacement_
    (
        coeffDict().getOrDefault<bool>
        ("separateTransRotPointDisplacement", false)
    ),
    //rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    scale_
    (
        IOobject
        (
            "motionScale_" + dict.dictName(),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    curTimeIndex_(-1),
    DoFs_(0)
{

    // Determining active degrees of freedom
    
    Vector<bool> linDirs(coeffDict().getOrDefault
    (
        "linDirs", Vector<bool>(1, 1, 1))
    );
    Vector<bool> rotDirs(coeffDict().getOrDefault
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

    // Adding linear degrees of freedom to DoFs_
    forAll(linDirs, i)
    {
        if (linDirs[i])
        {
            DoFs_.append(i);
        }
    }

    // Adding rotational degrees of freedom to DoFs_
    forAll(rotDirs, i)
    {
        if (rotDirs[i])
        {
            DoFs_.append(i+3);
        }
    }
   Info << "Active degrees of freedom: " << DoFs_ << endl;


   // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const pointMesh& pMesh = pointMesh::New(mesh);

//        pointPatchDist pDist(pMesh, patchSet_, points0_);
        pointPatchDist pDist(pMesh, patchSet_, (distStretch_ & points0_));

        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    (do_ - pDist.primitiveField())/(do_ - di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale_.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );

        pointConstraints::New(pMesh).constrain(scale_);
//        scale_.write();
    }

    // Reading settings for translational point displacement
    if (separateTransRotPointDisplacement_)
    {
        horDir1_ = coeffDict().get<vector>("horDir1");
        horDir2_ = coeffDict().get<vector>("horDir2");
        verDir_ = coeffDict().get<vector>("verDir");

        horPos1_ = coeffDict().get<scalarList>("horPos1");
        horPos2_ = coeffDict().get<scalarList>("horPos2");
        verPos_ = coeffDict().get<scalarList>("verPos");
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::floaterMotion&
Foam::floaterMeshMotionSolver::motion() const
{
    return motion_;
}


Foam::floaterMotion&
Foam::floaterMeshMotionSolver::motion()
{
    return motion_;
}


Foam::tmp<Foam::pointField>
Foam::floaterMeshMotionSolver::curPoints() const
{
    return points0_ + pointDisplacement_.primitiveField();
}


void Foam::floaterMeshMotionSolver::solve()
{

    if (mesh().nPoints() != points0_.size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0_.size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        motion_.storeOldState();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    // Body position and orientation updates were here but are now handled by
    // floaterMotion

    // Update the displacements
    updateDisplacement();


    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


void Foam::floaterMeshMotionSolver::updateDisplacement()
{

    point centreOfRotation(motion().centreOfRotation());

    if (separateTransRotPointDisplacement_)
    {
        centreOfRotation = motion().initialCentreOfRotation();
    }

    // Calculate the transformation septerion from the initial state
    septernion s
    (
        centreOfRotation - motion().initialCentreOfRotation(),
        quaternion(motion().Q().T() & motion().initialQ())
    );

    pointField& disp = pointDisplacement_.primitiveFieldRef();

    forAll(disp, pointi)
    {
        // Move non-stationary points
        if (scale_[pointi] > SMALL)
        {
            // Use solid-body motion where scale = 1
            if (scale_[pointi] > 1 - SMALL)
            {
                disp[pointi] += centreOfRotation 
                  + (
                        ((motion().Q()) & (motion().initialQ().T()))
                          &
                        (points0_[pointi] - motion().initialCentreOfRotation())
                    )
                  - points0_[pointi];
            }
            // Slerp septernion interpolation
            else
            {
                septernion ss(slerp(septernion::I, s, scale_[pointi]));

                disp[pointi] +=
                    motion().initialCentreOfRotation()
                  + ss.invTransformPoint
                    (
                        points0_[pointi]
                      - motion().initialCentreOfRotation()
                    )
                  - points0_[pointi];

                // Test of alternative displacement interpolation to check if
                // quaternions are to blame for non volume conservation.
                // Not the case.
                /*
                point newPos = centreOfRotation 
                  + (
                        ((motion().Q()) & (motion().initialQ().T()))
                          &
                        (points0_[pointi] - motion().initialCentreOfRotation())
                    );

                disp[pointi] = scale_[pointi]*newPos
                    + (1.0 - scale_[pointi])*points0_[pointi]
                    - points0_[pointi];
                */
            }
        }
    }

    if (separateTransRotPointDisplacement_)
    {
        const vector bodyDisp
        (
            motion().centreOfRotation() - motion().initialCentreOfRotation()
        );
        scalarField s(scale_.size(), 1.0);

        // Mesh point displacement due to x body displacement
        vector xhat(horDir1_/mag(horDir1_));
        calcTransDispScale(points0_, horPos1_, xhat, s);
        const scalar Lx(xhat & bodyDisp);
        disp += Lx * xhat * max
        (
            min(0.5 - 0.5*cos(Foam::constant::mathematical::pi*s), 1.0),
            0.0
        );

        // Mesh point displacement due to y body displacement
        s = 1.0;
        vector yhat(horDir2_/mag(horDir2_));
        calcTransDispScale(points0_, horPos2_, yhat, s);
        const scalar Ly(yhat & bodyDisp);
        disp += Ly * yhat * max
        (
            min(0.5 - 0.5*cos(Foam::constant::mathematical::pi*s), 1.0),
            0.0
        );

        // Mesh point displacement due to y body displacement
        s = 1.0;
        vector zhat(verDir_/mag(verDir_));
        calcTransDispScale(points0_, horPos1_, xhat, s);
        calcTransDispScale(points0_, horPos2_, yhat, s);
        calcTransDispScale(points0_, verPos_, zhat, s);
        const scalar Lz(zhat & bodyDisp);
        disp += Lz * zhat * max
        (
            min(0.5 - 0.5*cos(Foam::constant::mathematical::pi*s), 1.0),
            0.0
        );
    }
}


void Foam::floaterMeshMotionSolver::calcTransDispScale
(
    const pointField& X,
    const scalarList& pts,
    const vector& dir,
    scalarField& scale
) const
{
    const vector unitDir(dir/mag(dir));
    scale *= max(min(((X & unitDir) - pts[0])/(pts[1] - pts[0]), 1.0), 0.0);
    scale *= max(min((pts[3] - (X & unitDir)) / (pts[3] - pts[2]), 1.0), 0.0);
}


bool Foam::floaterMeshMotionSolver::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    IOdictionary dict
    (
        IOobject
        (
            "floaters",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            false
        )
    );

    // Create or get the sub-dictionary with the name from motion_.name()
    dictionary& floaterDict = dict.subDictOrAdd(motion_.name());
    motion_.write(floaterDict);

    return dict.regIOobject::write();
}


bool Foam::floaterMeshMotionSolver::read()
{
    if (pointDisplacementMotionSolver::read())
    {
        motion_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::labelList& Foam::floaterMeshMotionSolver::DoFs()
{
    return DoFs_;
}

// ************************************************************************* //
