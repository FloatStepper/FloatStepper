[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8146515.svg)](https://doi.org/10.5281/zenodo.8146515)

# FloatStepper 

OpenFOAM extension module for fluid-rigid body coupling.

Works for two-phase flows and also in single-phase mode by setting alpha field to 1 internally and on boundaries.

FloatStepper solves the added mass instability problem by calculating the added mass matrix, $A$, of the rigid body and exploiting $A$ in the updating of the body motion.

The basic idea is the following decomposition of the net force (and torque) on a body into an added mass term and everything else,

$\mathbf F = \mathbf F_\textrm{other} - A\mathbf a$, 

where $\mathbf a$ is the body acceleration, and $\mathbf F_\textrm{other}$ contains all forces not proportional to $\mathbf a$.
Inserting this into Newton's 2nd law, we have

$M\mathbf a = \mathbf F_\textrm{other} - A\mathbf a \Rightarrow \mathbf a = (M + A)^{-1} \mathbf F_\textrm{other}$

FloatStepper finds $\mathbf F_\textrm{other}$ by preceeding the actual time step by a test time step with zero body acceleration.\
The added mass matrix components are numerically calculated as the force (and torque) per (linear and angular) acceleration of the body.

A 15-minute presentation of FloatStepper from the 18th OpenFOAM Workshop can be found [here](https://youtu.be/Nn3Zl1jnr5U)

With **FloatStepper** you no longer need outer correctors or acceleration relaxation to stabilise coupling.

FloatStepper is not necessarily faster than the existing sixDoFRigidBodyMotion library in OpenFOAM, but should be more stable for light bodies in heavy fluids.

In terms of code maturity, it is at the proof-of-concept stage so users should expect (and report) bugs and also be willing to accept some syntax and API changes in future updates.

## Solver
*floatStepper*

## Library
*floaterMotion* (code structure based on *sixDoFRigidBodyMotion*)

## Requirements
OpenFOAM-v2206 to v2312. 
Only openfoam.com versions supported.

## Installation
1. Source OpenFOAM-v2XYZ
2. Execute `./Allwmake` from FloatStepper main directory
3. (optional) To install MoorDyn coupling go to thirdparty/MoorDyn and execute
   the Allwmake script found there.

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
    - `rotDirs (0 1 0);`
  
    which means translation only along x- and z-axes and rotation only around the y-axis (default is all DoF's active).
    The integer parameter *MaddUpdateFreq* (default 1) determines how often the (computationally expensive) added mass update is done, so setting *MaddUpdateFreq* e.g. to 3 means that the added mass matrix is only updated evert 3rd time step.
    sixDofRigidBodyMotionCoeffs subdicts *constraints{}*, *solver{}* and parameters *accelerationRelexation* and *accelerationDamping* will not be read or used by **FloatStepper**.
3.  Use *floaterVelocity* for U on floating object patches. Set `slip true;` to run with slip boundary condition (default is `no-slip`).
4.  Specify in `0.orig/uniform/floaterMotionState` the initial body position, orientation, linear and angular velocity and acceleration.

To test the MoorDyn installation, go to the run/moorDynBoxInWaves and execute the Allrun script. The mooring line parameters are located in the lines.txt file in run/moorDynBoxInWaves/Mooring. Mooring line forces will be written to files in the same folder.

## Author
Johan Roenby, STROMNING APS and Roskilde University

## Contributors
- Henrik Bredmose (conceptualisation)
- Sithik Aliyar (validation, verification, MoorDyn coupling)
- Henning Scheufler (code structure)

## Contributing
Please report bugs on the issue tracker of this repository or write to
johan[at]ruc[dot]dk

## Citing software
Roenby, J. (2024). FloatStepper (v0.0.3). Zenodo. https://doi.org/10.5281/zenodo.8146515

```bibtex
@software{roenby_2023_8146515,
  author       = {Roenby, Johan},
  title        = {FloatStepper},
  month        = jul,
  year         = 2023,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.8146515},
  url          = {https://doi.org/10.5281/zenodo.8146515}
}
```

## Publication
Roenby J, Aliyar S, Bredmose H., 2024 A robust algorithm for computational floating body dynamics. R. Soc. Open Sci. 11:231453. https://doi.org/10.1098/rsos.231453

```bibtex
@article{roenbyRobustAlgorithm2024,
	title = {A robust algorithm for computational floating body dynamics},
	volume = {11},
	issn = {2054-5703},
	url = {http://rsos.royalsocietypublishing.org/lookup/doi/10.1098/rsos.160405},
	doi = {10.1098/rsos.231453},
	language = {en},
	journal = {Royal Society Open Science},
	author = {Roenby, Johan and Aliyar, Sithik and Bredmose, Henrik},
	month = apr,
	year = {2024},
	pages = {231453}
}
```

## Known limitations
-   Currently only runs with morphing mesh (not overset mesh).
-   Currently only supports a single rigid body.
-   Has not been tested with turbulence modelling.
-   Currently only runs with isoAdvector (not MULES).
-   Plenty of field copying that can probably be avoided to increase efficiency.
-   Currently a simple time integration of the 6-DoF equations of motion are hardcoded in the *floaterMotion* class. The chosen method for Q update guarantees that it stays an orthogonal matrix to machine precision.

## Funding
Innovation Fund Denmark Grand Solution project FloatStep (8055-00075B).

Independen Research Fund Denmark, InterFlow project (9063-00018B).

Innovation Fund Denmark Grand Solution project FloatLab (2079-00042A).

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via www.openfoam.com, and owner of the 
OPENFOAM®  and OpenCFD®  trade marks.
