//Reset body state
forAll(bodySolvers, bsi)
{
    floaterMeshMotionSolver& bodySolver = bodySolvers[bsi];
    bodySolver.motion().setState(oldBodyStates[bsi]);
    //Note: would be nice to be able to do bodySolver.motion().setState(bodySolver.motion().oldState());
    //Info << "Body state after reset:" << endl;
    //bodySolver.motion().status();
}

// Reset mesh and swept face volumes
mesh.movePoints(mesh.oldPoints());
mesh.setPhi().ref() = mesh.phi().oldTime();
// Old mesh reset method (two steps back and one forward to recalculate swept
// face volumes):
//mesh.movePoints(oldOldPoints);
//mesh.movePoints(oldPoints);

alpha1 = alpha1.oldTime();
alpha2 = 1 - alpha1;
p_rgh = p_rgh.oldTime();
U = U.oldTime();
rho = rho.oldTime();
rhoPhi = rhoPhi.oldTime();
phi = phi.oldTime(); //Not needed if correctPhi
p = p.oldTime();
Uf() = Uf().oldTime();
