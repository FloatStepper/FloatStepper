Foam::myCorrectPhi
(
    U0,
    phi0,
    p_rgh0,
    surfaceScalarField("rAUf0", fvc::interpolate(rAU0)),
    geometricZeroField(),
    pimple
);

//#include "continuityErrs.H"
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