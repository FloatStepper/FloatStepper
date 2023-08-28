[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8146516.svg)](https://doi.org/10.5281/zenodo.8146516)


# FloatStepper 

Library for rigid body interaction with two phase flow.
Can also be used in single phase mode by setting alpha field to 1 internally and
on boundaries.
**FloatStepper** solves the added mass instability problem by calculating the added mass matrix of the rigid body and exploiting that in the updating of the body motion. 
With **FloatStepper** you no longer need outer correctors or acceleration relaxation to stabilise coupling. 

## Solver
*floatStepper*

## Library
*floaterMotion* (code structure based on *sixDoFRigidBodyMotion*)

## Requirements
OpenFOAM-v2206 (may compile with later versions as well).

## Installation
1. Source OpenFOAM-v2206
2. Execute `./Allwmake` from FloatStepper main directory

## Testing
Execute `./Allrun` in `run/circleFallingIntoWater`

Note that case syntax in slightly different from *sixDoFRigidBodyMotion* cases.
In particular:
1.  Set application to `floatStepper` in `system/controlDict`.
2.  Use *dynamicFloaterMotionSolversFvMesh* in `constant/dynamicMeshDict`.
    See examples of usage: 
    - [`run/circleFallingIntoWater/constant/dynamicMeshDict`](run/circleFallingIntoWater/constant/dynamicMeshDict).
    - [`run/mooredBoxInWaves/constant/dynamicMeshDict`](run/mooredBoxInWaves/constant/dynamicMeshDict).
    
    Restraints and their syntax is like *sixDoFRigidBodyMotion* but class names are changed, where "sixDoFRigidBody" is replaced with "floater".
    Constraints are specified with two vectors e.g.
    - `linDirs (1 0 1);`
    - `rotDirs (0 1 1);`
  
    which means translation only along x- and z-axes and rotation only around the y-axis (default is all DoF's active).
    The integer parameter *MaddUpdateFreq* (default 1) determines how often the (computationally expensive) added mass update is done, so setting *MaddUpdateFreq* e.g. to 3 means that the added mass matrix is only updated
    evert 3rd time step.
    Subdictionaries *constraints{}*, *solvers{}* and parameters *accelerationRelexation* and *accelerationDamping* will not be read or used by **FloatStepper**.
3.  Use *floaterVelocity* for U on floating object patches. Set `slip true;` to run with slip boundary condition (default is `no-slip`).
4.  Specify in `0.orig/uniform/floaterMotionState` the initial body position, orientation, linear and angular velocity and acceleration.

## Author
Johan Roenby, STROMNING APS and Roskilde University

## Contributors
- Henrik Bredmose (conceptualisation)
- Sithik Aliyar (validation, verification, moorDyn coupling (to appear))
- Henning Scheufler (code structure)

## Contributing
Please report bugs on the issue tracker of the repository or write to
johan@ruc.dk.

## Citing software
Roenby, Johan, (2023). FloatStepper-v0.0.1 [Computer Software]. Zenodo.
https://www.doi.org/10.5281/zenodo.8146516

```bibtex
@software{isoadvector_2023_8146516,
  author       = {Roenby, J.},
  title        = {FloatStepper},
  month        = jul,
  year         = 2023,
  publisher    = {Zenodo},
  version      = {v0.0.1},
  doi          = {10.5281/zenodo.8146516},
  url          = {https://doi.org/10.5281/zenodo.8146516}
}
```

## Known limitations
-   Currently only runs with morphing mesh (not overset mesh).
-   Currently only supports a single rigid body.
-   Has not been tested with turbulence modelling.
-   Currently only runs with isoAdvector (not MULES).
-   Plenty of field copying that can probably be avoided to increase efficiency.
-   Currently a simple time integration of the 6-DoF equations of motion are hardcoded in the *floaterMotion* class. The chosen method for Q update guarantees that it stays an orthogonal matrix to machine precission.
-   Possibly buggy behaviour when Center of Mass and Centre of Rotation do not coincide (needs testing).

## Funding
Innovation Fund Denmark Grand Solution project FloatStep (8055-00075B).

Independen Research Fund Denmark, InterFlow project (9063-00018B).

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via www.openfoam.com, and owner of the 
OPENFOAM®  and OpenCFD®  trade marks.