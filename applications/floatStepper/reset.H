//Reset body state
bodySolver.motion().setState(oldBodyState);
//Note: would be nice to be able to do bodySolver.motion().setState(bodySolver.motion().oldState());
//Info << "Body state after reset:" << endl;
//bodySolver.motion().status();


// Reset mesh and swept face volumes
mesh.movePoints(mesh.oldPoints());
mesh.setPhi().ref() = mesh.phi().oldTime();
// Old mesh reset method (two steps back and one forward to recalculate swept
// face volumes):
//mesh.movePoints(oldOldPoints);
//mesh.movePoints(oldPoints);

// Reset fields
//alpha1 = alpha1.oldTime(); //Reset does not work
//alpha2 = alpha2.oldTime(); //Reset does not work
p_rgh = p_rgh.oldTime();
U = U.oldTime();
rho = rho.oldTime();
rhoPhi = rhoPhi.oldTime();
phi = phi.oldTime(); //Not needed if correctPhi
p = p.oldTime();
Uf0 = Uf0.oldTime();
