{    
    if (mesh.moving())
    {
        // Temporarily making U relative to mesh motion
        U -= fvc::reconstruct(mesh.phi());

        // Workaround because isoAdvection::advect() saves wrong valuoes of 
        // alpha1.oldTime() when mesh is moving - part I
        alpha1.oldTime();
        alpha1.oldTime().primitiveFieldRef() *= (mesh.Vsc0()/mesh.Vsc());
    }

    // Updating alpha1
    #include "alphaSuSp.H"

    advector.advect(Sp, Su);

    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());

        // Workaround because isoAdvection::advect() saves wrong valuoes of
        // alpha1.oldTime() when mesh is moving - part II
        alpha1.oldTime().primitiveFieldRef() *= (mesh.Vsc()/mesh.Vsc0());
    }

    #include "rhofs.H"
    rhoPhi = advector.getRhoPhi(rho1f, rho2f);

    alpha2 = 1.0 - alpha1;
    mixture.correct();

}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
