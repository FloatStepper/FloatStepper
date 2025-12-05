// Note: For some reason alpha1 and alpha2 cannot be properly reset so we have
// to make copies of them for the zero acceleration time step.
volScalarField alpha10 = alpha1;
volScalarField alpha20 = alpha2;

// All other fields can be used for the 0-acceleration time step and then reset.
volVectorField& U0 = U;
surfaceScalarField& phi0 = phi;
isoAdvection advector0(alpha10, phi0, U0);
surfaceScalarField& rhoPhi0 = rhoPhi;
volScalarField& rho0 = rho;
volScalarField& p_rgh0 = p_rgh;
volScalarField& p0 = p;
surfaceVectorField& Uf0(Uf());//(Uf());
//volScalarField rAU0(rAU.ref());//(rAU());

// Accumulated displacement
// Copied and adapted from dynamicMotionSolverListMesh::update()
// Note: Only updates body motion solvers - not others.
mesh.update();
/*
{
    pointField disp(bodySolvers[0].newPoints() - mesh.points());

    for (label i = 1; i < bodySolvers.size(); i++)
    {
        disp += bodySolvers[i].newPoints() - mesh.points();
    }

    mesh.movePoints(mesh.points() + disp);

    volVectorField* Uptr = mesh.getObjectPtr<volVectorField>("U0");

    if (Uptr)
    {
        Uptr->correctBoundaryConditions();
    }
}
*/
/* //Implementation before multibody support
bodySolver.solve();
mesh.movePoints
(
    bodySolver.points0() + bodySolver.pointDisplacement().primitiveField()
);

U0.correctBoundaryConditions();
*/

// At this stage, the mesh points have been moved but all the fields are from
// the old mesh

uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        mesh.time().constant(),
        mesh,
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
volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

//MRF.update();

if (correctPhi)
{
    // Calculate absolute flux
    // from the mapped surface velocity
    phi0 = mesh.Sf() & Uf0;

    #include "correctPhi0.H"

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi0, U0);

    // mixture holds references to U an phi so to use mixture in 0-acceleration
    // time step we must temporarily overwrite these with U0 and phi0.
    // Note: mixture.correct is called again in alphaEqn0.H so needed here?
    mixture.correct();
}
/*
if (checkMeshCourantNo)
{
    #include "meshCourantNo.H"
}
*/