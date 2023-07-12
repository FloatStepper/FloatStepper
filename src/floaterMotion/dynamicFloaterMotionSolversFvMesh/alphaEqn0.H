/*---------------------------------------------------------------------------*\
    Module Name:       FloatStepper
    Description:       OpenFOAM extension module for fluid-rigid body coupling
    License:           GNU General Public License (GPL) version 3
    Copyright:         2023 Johan Roenby, STROMNING APS
\*---------------------------------------------------------------------------*/
{    
    if (fvMesh::moving())
    {
        U0 -= fvc::reconstruct(fvMesh::phi());
    }
    // Updating alpha0
    #include "alphaSuSp.H"
    advector0.advect(Sp, Su);

    if (fvMesh::moving())
    {
        U0 += fvc::reconstruct(fvMesh::phi());
    }

//    #include "rhofs.H"
    rhoPhi0 = advector0.getRhoPhi(rho1, rho2);

    alpha20 = 1.0 - alpha0;
//    mixture.correct();

}

Info<< "Phase-1 volume fraction = "
    << alpha0.weightedAverage(Vsc()).value()
    << "  Min(" << alpha0.name() << ") = " << min(alpha0).value()
    << "  Max(" << alpha0.name() << ") = " << max(alpha0).value()
    << endl;
