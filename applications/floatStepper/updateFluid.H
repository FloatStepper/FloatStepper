/*---------------------------------------------------------------------------*\
    Module Name:       FloatStepper
    Description:       OpenFOAM extension module for fluid-rigid body coupling
    License:           GNU General Public License (GPL) version 3
    Copyright:         2023 Johan Roenby, STROMNING APS
\*---------------------------------------------------------------------------*/
U.correctBoundaryConditions();

// --- Pressure-velocity PIMPLE corrector loop
while (pimple.loop())
{

    #include "alphaControls.H"
    #include "alphaEqnSubCycle.H"

    mixture.correct();

    #include "UEqn.H"

    // --- Pressure corrector loop
    while (pimple.correct())
    {
        #include "pEqn.H"
    }

    if (pimple.turbCorr())
    {
        turbulence->correct();
    }
}
