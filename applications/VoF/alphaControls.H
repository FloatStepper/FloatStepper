/*---------------------------------------------------------------------------*\
    Module Name:       FloatStepper
    Description:       OpenFOAM extension module for fluid-rigid body coupling
    License:           GNU General Public License (GPL) version 3
    Copyright:         2023 Johan Roenby, STROMNING APS
\*---------------------------------------------------------------------------*/
const dictionary& alphaControls = mesh.solverDict(alpha1.name());

//label nAlphaCorr(alphaControls.get<label>("nAlphaCorr"));

label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

/*
bool MULESCorr(alphaControls.lookupOrDefault("MULESCorr", false));

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed
bool alphaApplyPrevCorr
(
    alphaControls.lookupOrDefault("alphaApplyPrevCorr", false)
);

// Isotropic compression coefficient
scalar icAlpha
(
    alphaControls.lookupOrDefault<scalar>("icAlpha", 0)
);

// Shear compression coefficient
scalar scAlpha
(
    alphaControls.lookupOrDefault<scalar>("scAlpha", 0)
);
*/
