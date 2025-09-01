/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

//    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho0, U0) + fvm::div(rhoPhi0, U0, "div(rhoPhi,U)")
//      + MRF.DDt(rho0, U0)
      + turbulence->divDevRhoReff(rho0, U0) //NOTE: Uses alpha1, not alpha10 as it should. So term is lagged in time!!!
//     ==
//        fvOptions(rho, U)
    );

//    UEqn.relax();

//    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce() //NOTE: Uses alpha1, not alpha10 as it should. So term is lagged in time!!!
                  - ghf*fvc::snGrad(rho0)
                  - fvc::snGrad(p_rgh0)
                ) * mesh.magSf()
            )
        );

//        fvOptions.correct(U);
    }
