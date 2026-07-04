/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#include "exponentialMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exponentialMap, 0);
    addToRunTimeSelectionTable(floaterIntegrator, exponentialMap, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::exponentialMap::exponentialMap(const word& name, const dictionary& dict)
:
    floaterIntegrator(name, dict)
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

void Foam::exponentialMap::integrate
(
    floaterMotionState& motionState,
    const scalarField& dvwdt,
    const scalar deltaT
)
{
    // Store current state
    vector a0 = motionState.a();
    vector alpha0 = motionState.domegadt();
    vector a(dvwdt[0], dvwdt[1], dvwdt[2]);
    vector alpha(dvwdt[3], dvwdt[4], dvwdt[5]);
    motionState.a() = a;
    motionState.domegadt() = alpha;

    // Linear motion: update position with quadratic term (2nd-order kinematics)
    vector x0 = motionState.centreOfRotation();
    vector v0 = motionState.v();
    motionState.centreOfRotation() = x0 + v0*deltaT + 0.5*a0*deltaT*deltaT;
    motionState.v() = v0 + a0*deltaT;

    // Angular motion: exponential map with 2nd-order angular acceleration
    vector omega0 = motionState.omega();
    vector alpha_dt2 = 0.5*alpha0*deltaT*deltaT;

    // Effective rotation vector accounting for angular acceleration
    vector omega_eff = omega0*deltaT + alpha_dt2;
    scalar magw = mag(omega_eff);

    if (magw > SMALL)
    {
        // Rodrigues formula with effective rotation angle
        vector omega_hat = omega_eff / magw;
        tensor B
        (
            0, -omega_hat[2], omega_hat[1],
            omega_hat[2], 0, -omega_hat[0],
            -omega_hat[1], omega_hat[0], 0
        );

        tensor Qnew =
        (
            tensor::I + B*Foam::sin(magw) + (B & B)*(1 - Foam::cos(magw))
        ) & motionState.Q();
        motionState.Q() = Qnew;
    }

    // Update angular velocity
    motionState.omega() = omega0 + alpha0*deltaT;
}

// ************************************************************************* //
