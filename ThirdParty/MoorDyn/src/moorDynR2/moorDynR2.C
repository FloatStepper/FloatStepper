/*--------------------------------------------------- -----------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*--------------------------------------------------- -----------------------*/
/*---------------------------------------------------------------------------*\
References

    Chen, H., & Hall, M. (2022). CFD simulation of floating body motion with
        mooring dynamics: Coupling MoorDyn with OpenFOAM. Applied Ocean
        Research, 124, 103210. https://doi.org/10.1016/j.apor.2022.103210

\*---------------------------------------------------------------------------*/

#include "moorDynR2.H"
#include "addToRunTimeSelectionTable.H"
#include "rigidBodyMotion.H"
//#include "Time.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "error.H"
#include "quaternion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rigidBodyMotionRestraints
{
    defineTypeNameAndDebug(moorDynR2, 0);

    addToRunTimeSelectionTable
    (
        rigidBodyMotionRestraint,
        moorDynR2,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodyMotionRestraints::moorDynR2::moorDynR2
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    rigidBodyMotionRestraint(name, sDoFRBMRDict)
{
    read(sDoFRBMRDict);
    initialized_ = false;
    moordyn_backup_.t = 0.0;
    moordyn_backup_.data = nullptr;

    // Create the MoorDyn system
    if (Pstream::master())
    {
        int moordyn_err = MOORDYN_SUCCESS;
        moordyn_ = MoorDyn_Create("Mooring/lines_v2.txt");
        if (!moordyn_) {
            FatalError << "MoorDyn v2 cannot be created!" << exit(FatalError);
        }
        // In this release just a single floating body is accepted. In the
        // future more bodies can be added, providing info about the body to
        // be modelled
        unsigned int n;
        moordyn_err = MoorDyn_NCoupledDOF(moordyn_, &n);
        if ((moordyn_err != MOORDYN_SUCCESS) || (n != 6)) {
            FatalError << "6 Degrees Of Freedom were expected "
                    << "on the MoorDyn definition file" << exit(FatalError);
        }
        moordyn_err = MoorDyn_GetNumberBodies(moordyn_, &n);
        if ((moordyn_err != MOORDYN_SUCCESS) || !n) {
            FatalError << "At least one body was expected "
                    << "on the MoorDyn definition file" << exit(FatalError);
        }
        moordyn_body_ = NULL;
        for (unsigned int i = 0; i < n; i++) {
            MoorDynBody body = MoorDyn_GetBody(moordyn_, i + 1);
            if (!body) {
                FatalError << "Failure getting the MoorDyn body " << i + 1
                        << exit(FatalError);
            }
            int t;
            moordyn_err = MoorDyn_GetBodyType(body, &t);
            if (moordyn_err != MOORDYN_SUCCESS) {
                FatalError << "Failure geeting the body " << i + 1
                        << " type" << exit(FatalError);
            }
            if (t == -1) {
                // Coupled body, see Body.hpp:143
                moordyn_body_ = body;
                break;
            }
        }
        if (!moordyn_body_) {
            FatalError << "No coupled body could be found" << exit(FatalError);
        }
    }
    
    Info << "Create moorDynR2 using MoorDyn v2." << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rigidBodyMotionRestraints::moorDynR2::~moorDynR2()
{
    if (initialized_)
    {
        // Close MoorDyn call
        MoorDyn_Close(moordyn_);
        if (moordyn_backup_.data)
            free(moordyn_backup_.data);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::rigidBodyMotionRestraints::moorDynR2::restrain
(
    const rigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    int moordyn_err = MOORDYN_SUCCESS;

    scalar deltaT = motion.time().deltaTValue();
    scalar t = motion.time().value();
    scalar tprev = t - deltaT;

    point CoM = motion.centreOfMass();
    vector rotationAngle
    (
       quaternion(motion.orientation()).eulerAngles(quaternion::XYZ)
    );

    vector v = motion.v();
    vector omega = motion.omega();

    double X[6], XD[6];
    for (int ii=0; ii<3; ii++)
    {
       X[ii] = CoM[ii];
       X[ii+3] = rotationAngle[ii];
       XD[ii] = v[ii];
       XD[ii+3] = omega[ii];
    }

    if (!initialized_)
    {
        moordyn_err = MoorDyn_Init(moordyn_, X, XD);
        if (moordyn_err != MOORDYN_SUCCESS) {
            FatalError << "MoorDyn could not be initialized"
                       << exit(FatalError);
        }
        Info<< "MoorDyn module initialized!" << endl;
        initialized_ = true;
        save_mooring(tprev);
    } else if (tprev - moordyn_backup_.t >= 1.e-3 * deltaT) {
        // We have successfully advanced forward in time
        save_mooring(tprev);
        Info<< "MoorDyn module saved at t = " << tprev << " s" << endl;
    } else {
        // We are repeating the same time step because the implicit scheme
        load_mooring();
        Info<< "MoorDyn module restored to t = " << moordyn_backup_.t << " s" << endl;
    }

    double Flines[6] = {0.0};

    // Call LinesCalc() to obtain forces and moments, Flines(1x6)
    // LinesCalc(double X[], double XD[], double Flines[], double* t_in, double* dt_in)

    Info << "X[6]: " << vector(X[0], X[1], X[2]) << ", "
         << vector(X[3], X[4], X[5]) << endl;

    Info << "XD[6]: " << vector(XD[0], XD[1], XD[2]) << ", "
         << vector(XD[3], XD[4], XD[5]) << endl;

    moordyn_err = MoorDyn_Step(moordyn_, X, XD, Flines, &tprev, &deltaT);
    if (moordyn_err != MOORDYN_SUCCESS) {
        FatalError << "Error computing MoorDyn step " << tprev
              << "s -> " << tprev + deltaT << "s" << exit(FatalError);
    }
    
    // Info << "Mooring [force, moment] = [ "
    //      << Flines[0] << " " << Flines[1] << " " << Flines[2] << ", "
    //      << Flines[3] << " " << Flines[4] << " " << Flines[5] << " ]"
    //      << endl;

    for(int i=0;i<3;i++)
    {
        restraintForce[i] = Flines[i];
        restraintMoment[i] = Flines[i+3];
    }

    // Since moment is already calculated by LinesCalc, set to
    // centreOfRotation to be sure of no spurious moment
    restraintPosition = motion.centreOfRotation();

    //if (motion.report())
    {
        Info << t << ": force " << restraintForce
             << ", moment " << restraintMoment
             << endl;
    }
}


bool Foam::rigidBodyMotionRestraints::moorDynR2::read
(
    const dictionary& sDoFRBMRDict
)
{
    rigidBodyMotionRestraint::read(sDoFRBMRDict);

    return true;
}


void Foam::rigidBodyMotionRestraints::moorDynR2::write
(
    Ostream& os
) const
{
}


// ************************************************************************* //
