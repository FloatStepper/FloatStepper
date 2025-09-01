{    
    if (mesh.moving())
    {
        U0 -= fvc::reconstruct(mesh.phi());
    }
    // Updating alpha10
    #include "alphaSuSp.H"
    advector0.advect(Sp, Su);

    if (mesh.moving())
    {
        U0 += fvc::reconstruct(mesh.phi());
    }

//    #include "rhofs.H"
    rhoPhi0 = advector0.getRhoPhi(rho1, rho2);

    alpha20 = 1.0 - alpha10;
    // mixture holds references to U an phi so to use mixture in 0-acceleration
    // time step we must temporarily overwrite these with U0 and phi0.
/*
    volVectorField Utmp(U);
    surfaceScalarField phitmp(phi);
    U = U0;
    phi = phi0;
*/
    mixture.correct();
/*
    U = Utmp;
    phi = phitmp;
*/
}

Info<< "Phase-1 volume fraction = "
    << alpha10.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha10.name() << ") = " << min(alpha10).value()
    << "  Max(" << alpha10.name() << ") = " << max(alpha10).value()
    << endl;
