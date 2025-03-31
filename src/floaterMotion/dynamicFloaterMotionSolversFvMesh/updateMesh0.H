const immiscibleIncompressibleTwoPhaseMixture& mixture
    = lookupObject<immiscibleIncompressibleTwoPhaseMixture>("transportProperties");

const volScalarField& alpha1(mixture.alpha1());
const volScalarField& alpha2(mixture.alpha2());
//const volScalarField& alpha1 = lookupObject<volScalarField>("alpha.water");
//const volScalarField& alpha2 = lookupObject<volScalarField>("alpha.air");
const volScalarField& rho = lookupObject<volScalarField>("rho");
const surfaceScalarField& rhoPhi = lookupObject<surfaceScalarField>("rhoPhi");
const volVectorField& U = lookupObject<volVectorField>("U");
const surfaceScalarField& phif = lookupObject<surfaceScalarField>("phi");
const volScalarField& p_rgh = lookupObject<volScalarField>("p_rgh");
const volScalarField& p = lookupObject<volScalarField>("p");
const surfaceVectorField& Uf = lookupObject<surfaceVectorField>("Uf");
const volScalarField& rAU = lookupObject<volScalarField>("rAU");

volScalarField alpha0(alpha1);
volScalarField alpha20(alpha2);
volVectorField U0("U0", U);
surfaceScalarField phi0(phif);
isoAdvection advector0(alpha0, phi0, U0);
surfaceScalarField rhoPhi0(rhoPhi);
volScalarField rho0(rho);
volScalarField p_rgh0(p_rgh);
volScalarField p0("p0", p);
surfaceVectorField Uf0(Uf);//(Uf());
volScalarField rAU0(rAU);//(rAU());

bodySolver.solve();
disp += bodySolver.pointDisplacement().primitiveField();
fvMesh::movePoints(bodySolver.points0() + disp);
disp -= bodySolver.pointDisplacement().primitiveField();
//movePoints(bodySolver.newPoints());

U0.correctBoundaryConditions();

// At this stage, the mesh point have been moved but all the fields are from the
// old mesh

uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        time().constant(),
        *this,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimLength, Zero)
);

dimensionedVector g("g", dimAcceleration, Zero);
g = rho.time().lookupObject<uniformDimensionedVectorField>("g");

dimensionedScalar ghRef
(
    mag(g.value()) > SMALL
    ? g & (cmptMag(g.value())/mag(g.value()))*hRef
    : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
);
volScalarField gh("gh", (g & C()) - ghRef);
surfaceScalarField ghf("ghf", (g & Cf()) - ghRef);

//MRF.update();

pimpleControl& pimple =
    const_cast<pimpleControl&>
    (
        lookupObject<pimpleControl>("solutionControl")
    );

bool correctPhi
(
    pimple.dict().getOrDefault("correctPhi", dynamic())
);

if (correctPhi)
{
    // Calculate absolute flux
    // from the mapped surface velocity
    phi0 = Sf() & Uf0;

    #include "correctPhi0.H"

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi0, U0);

//    mixture.correct();
}
/*
if (checkMeshCourantNo)
{
    #include "meshCourantNo.H"
}
*/