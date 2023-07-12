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
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "floaterMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forceContributions.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"
#include "isoAdvection.H"
#include "CorrectPhi.H"
#include "pimpleControl.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(floaterMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        floaterMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotionSolver::floaterMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    motion_
    (
        coeffDict(), //Inherited from motionSolver - reads from dynamicMeshDict
        IOobject
        (
            "floaterMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "floaterMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : coeffDict(),
        mesh.time()
    ),
    patches_(coeffDict().get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(coeffDict().get<scalar>("innerDistance")),
    do_(coeffDict().get<scalar>("outerDistance")),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    scale_
    (
        IOobject
        (
            "motionScale",
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

        pointPatchDist pDist(pMesh, patchSet_, points0());

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
        scale_.write();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::floaterMotion&
Foam::floaterMotionSolver::motion() const
{
    return motion_;
}


Foam::floaterMotion&
Foam::floaterMotionSolver::motion()
{
    return motion_;
}


Foam::tmp<Foam::pointField>
Foam::floaterMotionSolver::curPoints() const
{
    return points0() + pointDisplacement_.primitiveField();
}


void Foam::floaterMotionSolver::solve()
{
//    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        motion_.storeOldState();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    // Removed body update - handled in floaterFoam


// if (false)
// {
// /*
//     // Grab handle to mesh motion solver
//     floaterMotionSolver& bodySolver =
//         const_cast<floaterMotionSolver&>
//         (
//             mesh().lookupObject<floaterMotionSolver>("dynamicMeshDict")
//         );
// */

//     const scalar deltaT = time().deltaTValue();

//     const fvMesh& mesh_ = dynamic_cast<const fvMesh&>(mesh());

//     if (time().timeIndex() == time().startTimeIndex() + 1)
//     {
//         Info << "First time step - not calculating added mass" << endl;
//         vector a = motion().state().a();
//         vector alpha = motion().state().domegadt();
//         scalarField dvwdt(6, 0);
//         dvwdt[0] = a[0], dvwdt[1] = a[1], dvwdt[2] = a[2];
//         dvwdt[3] = alpha[0], dvwdt[4] = alpha[1], dvwdt[5] = alpha[2];
//         motion().updateFloaterState(dvwdt, deltaT);
//     }
//     else
//     {
//         pointField oldPoints(mesh_.points());
//         pointField oldOldPoints(mesh_.oldPoints());
//         // oldMeshPhi used in alternative way to reset mesh motion (see below)
//         // surfaceScalarField oldMeshPhi(mesh.phi());
//         const floaterMotionState oldBodyState(motion().state());
//         const scalarField noAcceleration(6,0);
//         motion().setAcceleration(Zero);
//         motion().setAngularAcceleration(Zero);
//         motion().updateFloaterState(noAcceleration, deltaT);
//         Info << "Body state after zero acceleration time step:" << endl;
//         motion().status();

//         // Update mesh points in accordance with 0-acceleration body motion
//         fvMesh& nonConstPMesh = const_cast<fvMesh&>(mesh_);
//         #include "updateMesh0.H"
//         // Calculate fluid motion corresponding to 0-acceleration body motion
//         #include "updateFluid0.H" //alphaEqn, UEqn and Piso loop
//         // Recording fluid force associated with 0-acceleration motion
//         vector F0(Zero);
//         vector tau0(Zero);
//         motion().calcForceAndTorque(rho0, p0, U0, F0, tau0);
//         Info << "F0 = " << F0 << ", tau0 = " << tau0 << endl;

//         // Resetting mesh going first two steps back, then one forward to
//         // regain correct mesh.phi()

//         nonConstPMesh.movePoints(oldOldPoints);
//         nonConstPMesh.movePoints(oldPoints);
//         // Alternative approach with storage of old meshPhi
//         // mesh.setPhi() = oldMeshPhi;
//         // mesh.movePoints(oldPoints);

//         //Reset body state
//         motion().setState(oldBodyState);
//         Info << "Body state after reset:" << endl;
//         motion().status();

//         // Calculating the 6 columns of the added mass tensor
// //        const labelList DoFs = DoFs();

//         if (motion().MaddUpdateTime())
//         {
//             nonConstPMesh.moving(false);
//             motion().calcAddedMass(mesh_, DoFs());
//             nonConstPMesh.moving(true);
//         }
//         // Calculating 6 element acceleration vector 
//         scalarField dvwdt = motion().calcAcceleration(F0, tau0, DoFs());
//         // Updating body state to new time
//         motion().updateFloaterState(dvwdt, deltaT);
//         Info << "Body state in real time step:" << endl;
//         motion().status();

//     }
// }

    // Update the displacements
    pointDisplacement_.primitiveFieldRef() =
        motion_.transform(points0(), scale_) - points0();

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


bool Foam::floaterMotionSolver::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    IOdictionary dict
    (
        IOobject
        (
            "floaterMotionState",
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
}


bool Foam::floaterMotionSolver::read()
{
    if (displacementMotionSolver::read())
    {
        motion_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::labelList& Foam::floaterMotionSolver::DoFs()
{
    return DoFs_;
}

// ************************************************************************* //
