    Info<< "Creating velocity potential Pot\n" << endl;

    volScalarField Pot
    (
        IOobject
        (
            "Pot",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Pot", dimless, 0),
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

