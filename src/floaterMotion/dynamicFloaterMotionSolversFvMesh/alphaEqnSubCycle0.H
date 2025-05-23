if (pimple.nCorrPIMPLE() > 1)
{
    // If nOuterCorrectors > 1 then for all but the first loop the advection
    // of alpha is done using an average, 0.5*phi+0.5*phiNew where phi is
    // the flux at the beginning of the time step and phiNew is the flux
    // estimate at the end of the time step from the previous outer
    // iteration. Similarly we use 0.5*U + 0.5*UNew in later iterations.
    if (pimple.firstIter())
    {
        // To recalculate the alpha1 update in subsequent iterations, we
        // must store its current value before overwriting with the new
        // value
        alpha0.storePrevIter();
        // Storing initial phi and U for use in later outer iterations.
//        phi.storePrevIter();
//        U.storePrevIter();
//        Info << "store prevIter " << endl;
    }
    else
    {
        // Resetting alpha1 to value before advection in first PIMPLE
        // iteration.
//        Info << "set to prev Iter " << endl;
        alpha0 = alpha0.prevIter();

        // Setting U and phi with which to advect interface.
//        U = 0.5*U.prevIter() + 0.5*U;
//        phi = 0.5*phi.prevIter() + 0.5*phi;
    }
}

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = time().deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            time().timeName(),
            *this
        ),
        *this,
        dimensionedScalar(rhoPhi0.dimensions(), Zero)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha0, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn0.H"
        rhoPhiSum += (time().deltaT()/totalDeltaT)*rhoPhi0;
    }

    rhoPhi0 = rhoPhiSum;
}
else
{
    #include "alphaEqn0.H"
}

rho0 == alpha0*rho1 + alpha20*rho2;

if (!pimple.firstIter())
{
    // Resetting U and phi to value at latest iteration.
//    U = 2.0*U - U.prevIter();
//    phi = 2.0*phi - phi.prevIter();
}

