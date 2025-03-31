{
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    Foam::setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

/*
    if (correctPhi)
    {
        rAU0.ref() = 1.0/UEqn.A();
    }
    else
*/
    {
        rAU0 = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU0));
    volVectorField HbyA(constrainHbyA(rAU0*UEqn.H(), U0, p_rgh0));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
//      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
// If mesh.dynamic() ddtCorr(U,phi,Uf) redirects to ddtCorr(U, Uf())
      + fvc::interpolate(rho0*rAU0)*fvc::ddtCorr(U0, Uf0)
    );
//    MRF.makeRelative(phiHbyA);

    if (p_rgh0.needReference())
    {
        fvc::makeRelative(phiHbyA, U0);
        myAdjustPhi(phiHbyA, U0, p_rgh0);
        fvc::makeAbsolute(phiHbyA, U0);
    }

    surfaceScalarField phig
    (
        (
//            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho0)
        )*rAUf*magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
//    constrainPressure(p_rgh0, U0, phiHbyA, rAUf, MRF);
    constrainPressure(p_rgh0, U0, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh0) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(this->solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi0 = phiHbyA - p_rghEqn.flux();

            p_rgh0.relax();

            U0 = HbyA + rAU0*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U0.correctBoundaryConditions();
//            fvOptions.correct(U);
        }
    }

//    #include "continuityErrs.H"
/*
    {
        volScalarField contErr(fvc::div(phi0));

        scalar sumLocalContErr = runTime.deltaTValue()*
            mag(contErr)().weightedAverage(mesh.V()).value();

        scalar globalContErr = runTime.deltaTValue()*
            contErr.weightedAverage(mesh.V()).value();
        cumulativeContErr += globalContErr;

        Info<< "time step continuity errors : sum local = " << sumLocalContErr
            << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr
            << endl;
    }
*/
    // Correct Uf if the mesh is moving
//    fvc::correctUf(Uf0, U0, phi0);
    if (dynamic())
    {
        Uf0 = fvc::interpolate(U0);
        surfaceVectorField n(Sf()/magSf());
        Uf0 += n*(phi0/magSf() - (n & Uf0));
    }



    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi0, U0);

    p0 == p_rgh0 + rho0*gh;

    if (p_rgh.needReference())
    {
        p0 += dimensionedScalar
        (
            "p0",
            p0.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh0 = p0 - rho0*gh;
    }

    if (!correctPhi)
    {
//        rAU0.clear();
    }
}
