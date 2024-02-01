This directory contains base cases, case generation scripts (.sh) and figure
generation scripts (.py) to generate Figure 1-3 in "A robust algorithm for
computational floating body dynamics".

The case is a rigid disc in a large, circular domain, so the boundaries have
very small effect on the body motion. The fluid is a single phase ideal fluid,
i.e. no fluid interface and no viscosity. Gravity is set to -1 m/s2 along the
y-axis.

There is both a baseCase for the original sixDoFRigidBodyMotion body solver
using interIsoFoam, and for the floatStepper body solver.

Note that interIsoFoam supports two phase flows with a sharp fluid interface,
but can also be run in single phase mode by setting the volume fraction field
to 1. In this case the solver should behave as pimpleFoam.

Johan Roenby, October 2023