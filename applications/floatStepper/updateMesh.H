mesh.update();

Info << "Taking actual time step for fluid..." << endl;

gh = (g & mesh.C()) - ghRef;
ghf = (g & mesh.Cf()) - ghRef;

MRF.update();

if (correctPhi)
{
    // Calculate absolute flux
    // from the mapped surface velocity
    phi = mesh.Sf() & Uf();

    #include "correctPhi.H"

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi, U);

    mixture.correct();
}

if (checkMeshCourantNo)
{
    #include "meshCourantNo.H"
}
