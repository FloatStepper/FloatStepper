    Info<< "Creating velocity potential Phi\n" << endl;

    volScalarField Phi
    (
        IOobject
        (
            "Phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Phi", dimless, 0),
        "fixedGradient"
    );

    Info<< "Reading velocity field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

