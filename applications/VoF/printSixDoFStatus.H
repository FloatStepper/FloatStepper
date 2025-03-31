rigidBodyMotionSolver& bodySolver =
    const_cast<rigidBodyMotionSolver&>
    (
        mesh.lookupObject<rigidBodyMotionSolver>("dynamicMeshDict")
    );

// Variables needed for body and fluid momentum calculation
const rigidBodyMotion& bodyMotion(bodySolver.motion());
const rigidBodyMotionState& bodyState(bodyMotion.state());

bodyMotion.status();

/*
const vector& vBody = bodyState.v();
const scalar& mBody = bodyMotion.mass();
const scalarField& CV = mesh.V();
const vectorField& C = mesh.C();

vector PBody(vector::zero);
vector PFluid(vector::zero);
scalar EBody(0);
scalar EFluid(0);
scalar EpotFluid(0);
scalar EpotBody(0);


PBody = mBody*vBody;
PFluid = gSum(rho*CV*U);
vector Omega = bodyState.Q().T() & bodyState.omega();
EBody = 0.5*mBody*magSqr(vBody)
    + 0.5*(Omega & (bodyMotion.momentOfInertia() & Omega));
EFluid = 0.5*gSum(rho*CV*magSqr(U));
dimensionedVector gVec("gVec", dimAcceleration, Zero);
gVec = runTime.lookupObject<uniformDimensionedVectorField>("g");
vector zhat = -gVec.value();
zhat = zhat/(mag(zhat) + SMALL);
EpotFluid = gSum(rho*CV*mag(C & gVec.value()));
EpotBody = mBody*mag(bodyMotion.centreOfMass() & gVec.value());

Info << "Body linear momentum   : " << PBody << nl
    << "Fluid linear momentum  : " << PFluid << nl
    << "Total linear momentum  : " << PBody + PFluid << nl
    << "Body kinetic energy    : " << EBody << nl
    << "Fluid kinetic energy   : " << EFluid << nl
    << "Total kinetic energy   : " << EBody + EFluid << nl
    << "Fluid potential energy : " << EpotFluid << nl
    << "Body potential energy  : " << EpotBody << nl
    << "Total potential energy : " << EpotFluid + EpotBody << nl
    << "Total energy           : " << EpotFluid + EpotBody + EFluid + EBody
    << endl;
*/
