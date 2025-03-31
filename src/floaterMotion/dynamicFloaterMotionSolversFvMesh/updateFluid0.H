/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/
// --- Pressure-velocity PIMPLE corrector loop
while (pimple.loop())
{

//    #include "alphaControls.H"
    const dictionary& alphaControls = solverDict(alpha0.name());

    label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

    #include "alphaEqnSubCycle0.H"

//    mixture.correct();

    #include "UEqn0.H"

    // --- Pressure corrector loop
    while (pimple.correct())
    {
        #include "pEqn0.H"
    }

/*
    if (pimple.turbCorr())
    {
        turbulence->correct();
    }
*/
}
