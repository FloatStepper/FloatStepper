# FloatStepperPaperCases

This repo contains all OpenFOAM cases used to generate figures in the
FloatStepper method paper entitled "A robust algorithm for computational
floating body dynamics", published in Royal Society Open Science:

https://doi.org/10.1098/rsos.231453

All cases should contain a baseCase, and each baseCase should
contain the following:
- A complete Allrun script for meshing, running and postprocessing the case.
- An Allclean script cleaning all results, bringing the case folder back to 
  its exact original state.

The repository also contains bash scripts to generate and run the cases in the
paper from the base cases, and python scripts to generate the figures in the
paper from the cases.

To run the cases in the repo, the user must first install OpenFOAM-v2206 (or
later) from

https://www.openfoam.com/current-release

and install the floatStepper OpenFOAM extension from

https://github.com/FloatStepper/FloatStepper

Please report bugs using the issue tracker of the relevant repository.

Johan Roenby, October 2023
