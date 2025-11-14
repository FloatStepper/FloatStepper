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
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "floaterMotion.H"
#include "septernion.H"
#include "MatrixTools.H" //Only needed to printMatrix
#include "fvCFD.H" //Needed by calcAddedMass()
#include "forceContributions.H"
#include "uniformDimensionedFields.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(floaterMotion, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotion::floaterMotion(const Time& time)
:
    name_(),
    time_(time),
    motionState_(),
    motionState0_(),
    restraints_(),
//    constraints_(),
//    tConstraints_(tensor::I),
//    rConstraints_(tensor::I),
    initialCentreOfMass_(Zero),
    initialCentreOfRotation_(Zero),
    initialQ_(I),
    mass_(VSMALL),
    momentOfInertia_(symmTensor::one*VSMALL),
    momentOfInertiaRefPoint_(Zero),
    momentOfInertiaAxes_(Zero),
    Madd_({6, 6}, 0),
    MaddUpdateFreq_(1),
    F0_(Zero),
    tau0_(Zero),
    restraintForce_(Zero),
    restraintTorque_(Zero),
    patches_()
{}


Foam::floaterMotion::floaterMotion
(
    const dictionary& dict,
    const dictionary& stateDict,
    const Time& time
)
:
    name_(stateDict.dictName()),
    time_(time),
    motionState_(stateDict),
    motionState0_(),
    restraints_(),
//    constraints_(),
//    tConstraints_(tensor::I),
//    rConstraints_(tensor::I),
    initialCentreOfMass_(stateDict.get<point>("initialCentreOfMass")),
    initialCentreOfRotation_(stateDict.get<point>("initialCentreOfRotation")),
    initialQ_(stateDict.get<tensor>("initialOrientation")),
    mass_(stateDict.get<scalar>("mass")),
    momentOfInertia_(stateDict.get<symmTensor>("momentOfInertia")),
    momentOfInertiaRefPoint_(stateDict.getOrDefault<vector>("momentOfInertiaRefPoint", centreOfRotation())),
    momentOfInertiaAxes_(stateDict.getOrDefault<tensor>("momentOfInertiaAxes", orientation())),
    Madd_(stateDict.getOrDefault("Madd", SquareMatrix<scalar>({6, 6}, 0))),
    MaddUpdateFreq_(dict.getOrDefault("MaddUpdateFreq", 1)),
    F0_(stateDict.getOrDefault("F0", vector::zero)),
    tau0_(stateDict.getOrDefault("tau0", vector::zero)),
    restraintForce_(Zero),
    restraintTorque_(Zero),
    patches_(dict.get<wordRes>("patches"))
{

    Info << endl;
    Info << "Inertia for body: " << name_ << endl;
    Info << "mass: " << mass_ << endl;
    Info << "centreOfMass: " << centreOfMass() << endl;
    Info << "momentOfInertia: " << momentOfInertia_
        << " with respect to" << endl;
    Info << "momentOfInertiaRefPoint: " << momentOfInertiaRefPoint_ << " and "
        << endl;
    Info << "momentOfInertiaAxes: " << momentOfInertiaAxes_ << endl;

    // Transform momentOfInertia tensor to lab axes
    tensor J = (momentOfInertiaAxes_ & momentOfInertia_) & momentOfInertiaAxes_;
 
    // Calculate momentOfInertia with respect to centreOfMass
    vector d = momentOfInertiaRefPoint_ - centreOfMass();
    tensor JCofM = J - mass_*((d & d)*I - d*d);
 
    // Calculating momentOfInertia with respect to centreOfRotation
    d = centreOfRotation() - centreOfMass();
    tensor ICofR = JCofM + mass_*((d & d)*I - d*d);

    // Transforming momentOfInertia to body axes
    ICofR = (orientation().T() & ICofR) & orientation();

    // Changing momentOfInertia_ to be with respect to centre of rotation and in
    // body axes.
    momentOfInertia_ = symmTensor
    (
        ICofR.xx(), ICofR.xy(), ICofR.xz(),
        ICofR.yy(), ICofR.yz(), ICofR.zz()
    );
    momentOfInertiaRefPoint_ = centreOfRotation();
    momentOfInertiaAxes_ = orientation();

    Info << "momentOfInertia with respect to centre of rotation in body axes:"
        << endl;
    Info << momentOfInertia_ << endl;
    Info << endl;

    addRestraints(dict);

    // Save the old-time motion state
    motionState0_ = motionState_;
}


Foam::floaterMotion::floaterMotion
(
    const floaterMotion& rbm
)
:
    name_(rbm.name_),
    time_(rbm.time_),
    motionState_(rbm.motionState_),
    motionState0_(rbm.motionState0_),
    restraints_(rbm.restraints_),
//    constraints_(rbm.constraints_),
//    tConstraints_(rbm.tConstraints_),
//    rConstraints_(rbm.rConstraints_),
    initialCentreOfMass_(rbm.initialCentreOfMass_),
    initialCentreOfRotation_(rbm.initialCentreOfRotation_),
    initialQ_(rbm.initialQ_),
    mass_(rbm.mass_),
    momentOfInertia_(rbm.momentOfInertia_),
    Madd_(rbm.Madd_),
    MaddUpdateFreq_(rbm.MaddUpdateFreq_),
    F0_(rbm.F0_),
    tau0_(rbm.tau0_),
    restraintForce_(rbm.restraintForce_),
    restraintTorque_(rbm.restraintTorque_),
    patches_(rbm.patches_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//Foam::floaterMotion::~floaterMotion()
//{} // Define here (incomplete type in header)


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::floaterMotion::calcAddedMass
(
    const fvMesh& mesh,
    const labelList& DoFs
)
{
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );
    const dictionary& solvertDict = dynamicMeshDict.subDict("solvers");

    word bodyName;
    word motionSolverType;
    for (const entry& dEntry : solvertDict)
    {
        // Note! If more than one motionSolver, this ends up with the last one
        if (dEntry.isDict())
        {
            bodyName = dEntry.dict().dictName();
            const dictionary& bodyDict = solvertDict.subDict(bodyName);
            motionSolverType = bodyDict.get<word>("motionSolver");
        }
    }
    const dictionary& solverDict =
        solvertDict.subDict(bodyName).subDict(motionSolverType + "Coeffs");

    autoPtr<addedMass> addedM = addedMass::New
    (
        solverDict.get<word>("addedMassModel"),
        mesh,
        dynamicMeshDict
    );

    FixedList<bool,6> enabledDoFs(false);
    forAll (DoFs, DoFi)
    {
        enabledDoFs[DoFs[DoFi]] = true;
    }

    Madd_ = addedM->computeAddedMass
    (
        enabledDoFs,
        centreOfRotation(),
        Q(),
        patches_
    );

    // Changing to body frame
    Madd_ = changeFrame(Madd_, Q());
}


void Foam::floaterMotion::calcForceAndTorque
(
    const volScalarField& rho,
    const volScalarField& p,
    const volVectorField& U,
    vector& force,
    vector& torque
)
{

    dictionary forcesDict;

    forcesDict.add("type", functionObjects::forces::typeName);
    forcesDict.add("patches", patches_);
    forcesDict.add("rhoInf", 1.0);
    forcesDict.add("rho", rho.name());
    forcesDict.add("p", p.name());
    forcesDict.add("U", U.name());
    forcesDict.add("CofR", centreOfRotation());

    functionObjects::forceContributions f("forces", rho.time(), forcesDict);

    f.calcForcesMoments();

    dimensionedVector g("g", dimAcceleration, Zero);

    g = rho.time().lookupObject<uniformDimensionedVectorField>("g");

    applyRestraints();

    force = f.forceEff() + mass()*g.value() + restraintForce();

    torque = f.momentEff() + mass()*(momentArm() ^ g.value()) 
        + restraintTorque();

    F0_ = force;
    tau0_ = torque;
}


void Foam::floaterMotion::updateFloaterState
(
    const scalarField& dvwdt, 
    const scalar deltaT
)
{
/*
    // Preparation for runTime selection of 6DoF integrator
    autoPtr<RBMIntegrator> RBMUpdater = RBMIntegrator::New
    (
        dynamicMeshDict.get<word>("RBMIntegratorModel"),
        mesh,
        dynamicMeshDict
    );

    vector CoR0 = motionState_.CoR();
    vector Q0 = motionState_.Q();
    vector v0 = motionState_.v();
    vector omega0 = motionState_.omega();
    vector a0 = motionState_.a();
    vector alpha0 = motionState_.domegadt();

    vector a(dvwdt[0], dvwdt[1], dvwdt[2]);
    vector alpha(dvwdt[3], dvwdt[4], dvwdt[5]);

    RBMUpdater->integrateRBM
    (
        motionState_,
        dvwdt,
        deltaT
    );
*/

    // Updating body velocity and position
    vector a0 = motionState_.a();
    vector alpha0 = motionState_.domegadt();
    vector a(dvwdt[0], dvwdt[1], dvwdt[2]);
    vector alpha(dvwdt[3], dvwdt[4], dvwdt[5]);
    motionState_.a() = a;
    motionState_.domegadt() = alpha;

    label nSteps = 1e3;
    scalar dt = deltaT/nSteps;

    for (label n = 1; n <= nSteps; n++)
    {
        // Updating position according to Krysl & Endres 2005:
        vector x0 = motionState_.centreOfRotation();
        vector v0 = motionState_.v();
//        motionState_.centreOfRotation() = x0 + dt*(v0 + dt*0.5*a0);
        motionState_.centreOfRotation() = x0 + dt*v0;

        //Updating velocity
//        motionState_.v() = v0 + 0.5*dt*(a0 + a);
        motionState_.v() = v0 + dt*a0;
//        motionState_.v() = v0 + dt*a;

        // Updating body angular velocity and orientation
        vector omega0 = motionState_.omega();
//        motionState_.omega() = omega0 + 0.5*dt*(alpha0 + alpha);
        motionState_.omega() = omega0 + dt*alpha0;
//        motionState_.omega() = omega0 + dt*alpha;

        vector omega = omega0;
//        vector omega = omega0 + 0.5*dt*alpha0;
//        vector omega = 0.5*(omega0 + wnew);
//        scalar magw = mag(wnew);
        scalar magw = mag(omega);
        tensor B
        (
            0, -omega[2], omega[1],
            omega[2], 0, -omega[0],
            -omega[1], omega[0], 0
        );

        B /= (magw + SMALL);
        // Rodrigues rotation formula
        tensor Qnew =
        (
            tensor::I + B*Foam::sin(magw*dt) + (B & B)*(1 - Foam::cos(magw*dt))
        ) & motionState_.Q();
        motionState_.Q() = Qnew;
    }
}


Foam::scalarField Foam::floaterMotion::calcAcceleration
(
    const vector& F0,
    const vector& tau0,
    const labelList& DoFs
)
{
    // Body mass and inertia to effective mass matrix
    typedef SquareMatrix<scalar> SMatrix;
    SMatrix Mbody({6, 6}, 0);
    tensor Q = motionState_.Q();
    tensor I_lf = (Q & momentOfInertia_) & Q.T(); // lab frame inertia wrt CofR.
//    tensor Mtensor = mass()*(tensor::I); // Q mI QT = m Q Qt = mI
    vector marm = mass()*momentArm();
    tensor mXcm
    (
        0, -marm[2], marm[1],
        marm[2], 0, -marm[0],
        -marm[1], marm[0], 0
    );

//    Info << "mXcm = " << mXcm << ", momentArm = " << momentArm() << endl;
    Mbody.subMatrix(0,0,3,3) = mass()*(tensor::I); // Q mI QT = m Q Qt = mI
    Mbody.subMatrix(0,3,3,3) = mXcm.T();
    Mbody.subMatrix(3,0,3,3) = mXcm;
    Mbody.subMatrix(3,3,3,3) = I_lf;
//    Info << "Body mass matrix in lab frame:" << endl;
//    MatrixTools::printMatrix(Info, Mbody) << nl;

    // Defining RHS in the equation Meff*ddt[v, w] = [F0 - wxmXcmw, tau0 - wxIw]
    SMatrix Meff({6, 6}, 0);
    Meff = changeFrame(Madd_, Q.T()) + Mbody;
    scalarField RHS(6), dvwdt(6);
    vector omega = motionState_.omega();
    vector wxIw = omega ^ (I_lf & omega);
    vector wxmXcmw = omega ^ (mXcm & omega);
    RHS[0] = F0[0] - wxmXcmw[0];
    RHS[1] = F0[1] - wxmXcmw[1];
    RHS[2] = F0[2] - wxmXcmw[2];
    RHS[3] = tau0[0] - wxIw[0];
    RHS[4] = tau0[1] - wxIw[1];
    RHS[5] = tau0[2] - wxIw[2];
    dvwdt = RHS;

    // Clean inactive DoF
    const label nDoFs = DoFs.size();
    SMatrix Meff_reduced({nDoFs, nDoFs}, 0);
    scalarField dvwdt_reduced(nDoFs, 0);
    forAll(DoFs, di)
    {
        forAll(DoFs, dj)
        {
            Meff_reduced(di, dj) = Meff(DoFs[di], DoFs[dj]);
        }
        dvwdt_reduced[di] = dvwdt[DoFs[di]];
    }

    if (debug)
    {
        Info << "Effective mass matrix in lab frame:" << endl;
        MatrixTools::printMatrix(Info, Meff) << nl;

        Info << "Reduced effective mass matrix in lab frame:" << endl;
        MatrixTools::printMatrix(Info, Meff_reduced) << nl;

        Info << "Before LUsolve: dvwdt_reduced = " << dvwdt_reduced << endl;
    }

    // Solving Meff*ddt[v, w] = [F0 - wxmXcmw, tau0 - wxIw] with LU decomposition
    LUsolve(Meff_reduced, dvwdt_reduced);

    // Populating 6 DoF accelereation vector to be returned
    dvwdt = 0;
    forAll(dvwdt_reduced, di)
    {
        dvwdt[DoFs[di]] = dvwdt_reduced[di];
    }

    //Printing out Meff in body frame
    if (debug)
    {
        SMatrix Meff_bf = changeFrame(Meff, Q);
        Info << "After LUsolve: dvwdt_reduced = " << dvwdt_reduced << endl;
        Info << "Meff in body frame:" << endl;
        MatrixTools::printMatrix(Info, Meff_bf) << nl;
    }

    return dvwdt;
}

Foam::scalarSquareMatrix Foam::floaterMotion::changeFrame
(
    const scalarSquareMatrix& M,
    const tensor& Q
) const
{
    scalarSquareMatrix Qmat({6, 6}, 0);
    Qmat.subMatrix(0,0,3,3) = Q;
    Qmat.subMatrix(3,3,3,3) = Q;

    return (Qmat.T() * (M * Qmat));
}

bool Foam::floaterMotion::MaddUpdateTime() const
{
    return ((time_.timeIndex()-1) % MaddUpdateFreq_) == 0;
}

void Foam::floaterMotion::addRestraints
(
    const dictionary& dict
)
{
    if (dict.found("restraints"))
    {
        const dictionary& restraintDict = dict.subDict("restraints");

        label i = 0;

        restraints_.setSize(restraintDict.size());

        for (const entry& dEntry : restraintDict)
        {
            if (dEntry.isDict())
            {
                restraints_.set
                (
                    i++,
                    floaterMotionRestraint::New
                    (
                        dEntry.keyword(),
                        dEntry.dict()
                    )
                );
            }
        }

        restraints_.setSize(i);
    }
}


void Foam::floaterMotion::applyRestraints()
{
    if (restraints_.empty())
    {
        return;
    }

// Outcommented line below because it caused different subdomains to move 
// differently so rigid body changed shape. 
//    if (Pstream::master())
    {
        restraintForce_ = Zero;
        restraintTorque_ = Zero;

        forAll(restraints_, rI)
        {
            // Restraint position
            point rP = Zero;

            // Restraint force
            vector rF = Zero;

            // Restraint moment
            vector rM = Zero;

            // Accumulate the restraints
            restraints_[rI].restrain(*this, rP, rF, rM);

            // Update the force
            restraintForce_ += rF;

            // Moments are returned in global axes
            restraintTorque_ += rM + ((rP - centreOfRotation()) ^ rF);
        }
    }
}

/*
void Foam::floaterMotion::addConstraints
(
    const dictionary& dict
)
{
    if (dict.found("constraints"))
    {
        const dictionary& constraintDict = dict.subDict("constraints");

        label i = 0;

        constraints_.setSize(constraintDict.size());

        pointConstraint pct;
        pointConstraint pcr;

        for (const entry& dEntry : constraintDict)
        {
            if (dEntry.isDict())
            {
                constraints_.set
                (
                    i,
                    floaterMotionConstraint::New
                    (
                        dEntry.keyword(),
                        dEntry.dict(),
                        *this
                    )
                );

                constraints_[i].setCentreOfRotation(initialCentreOfRotation_);
                constraints_[i].constrainTranslation(pct);
                constraints_[i].constrainRotation(pcr);

                i++;
            }
        }

        constraints_.setSize(i);

        tConstraints_ = pct.constraintTransformation();
        rConstraints_ = pcr.constraintTransformation();

        Info<< "Translational constraint tensor " << tConstraints_ << nl
            << "Rotational constraint tensor " << rConstraints_ << endl;
    }
}
*/

void Foam::floaterMotion::status() const
{
    // Writing summary of body update
    Info << "--------------Body update--------------------" << endl;
    Info << "Centre of rotation  : " << motionState_.centreOfRotation() << endl;
    Info << "Velocity            : " << motionState_.v() << endl;
    Info << "Acceleration        : " << motionState_.a() << endl;
    Info << "Orientation         : " << motionState_.Q() << endl;
    Info << "Angular velocity    : " << motionState_.omega() << endl;
    Info << "Angular acceleration: " << motionState_.domegadt() << endl;
    Info << "---------------------------------------------" << endl;
}


Foam::tmp<Foam::pointField> Foam::floaterMotion::transform
(
    const pointField& initialPoints
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoints - initialCentreOfRotation_))
    );
}


// ************************************************************************* //
