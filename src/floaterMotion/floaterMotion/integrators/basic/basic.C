/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#include "basic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basic, 0);
    addToRunTimeSelectionTable(floaterIntegrator, basic, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::basic::basic(const word& name, const dictionary& dict)
:
    floaterIntegrator(name, dict)
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

void Foam::basic::integrate
(
    floaterMotionState& motionState,
    const scalarField& dvwdt,
    const scalar deltaT
)
{
    // Updating body velocity and position
    vector a0 = motionState.a();
    vector alpha0 = motionState.domegadt();
    vector a(dvwdt[0], dvwdt[1], dvwdt[2]);
    vector alpha(dvwdt[3], dvwdt[4], dvwdt[5]);
    motionState.a() = a;
    motionState.domegadt() = alpha;

    label nSteps = 1e3;
    scalar dt = deltaT/nSteps;

    for (label n = 1; n <= nSteps; n++)
    {
        // Updating position according to Krysl & Endres 2005:
        vector x0 = motionState.centreOfRotation();
        vector v0 = motionState.v();
        motionState.centreOfRotation() = x0 + dt*v0;

        // Updating velocity
        motionState.v() = v0 + dt*a0;

        // Updating body angular velocity and orientation
        vector omega0 = motionState.omega();
        motionState.omega() = omega0 + dt*alpha0;

        vector omega = omega0;
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
        ) & motionState.Q();
        motionState.Q() = Qnew;
    }
}

// ************************************************************************* //
