import fileinput, math, re

R1 = 1.0
R2 = 1.005
R3 = 50.0
z1 = -0.5
z2 = 0.5
nx = 200
ny = 60
nz = 1
grading = 60
Eps = 1.0e-3

R1OverSqrt2 = R1 / math.sqrt(2.0)
R2OverSqrt2 = R2 / math.sqrt(2.0)
R3OverSqrt2 = R3 / math.sqrt(2.0)

z1PlusEps = z1 + Eps
z1MinusEps = z1 - Eps
z2PlusEps = z2 + Eps
z2MinusEps = z2 - Eps

R1PlusEps = R1 + Eps
mR1PlusEps = - R1PlusEps

with open('system/meshParms', 'w') as meshParmFile:
    meshParmFile.write('R1    %.12f;\n' % R1);
    meshParmFile.write('R2    %.12f;\n' % R2);
    meshParmFile.write('R3    %.12f;\n' % R3);
    meshParmFile.write('R1OverSqrt2    %.12f;\n' % R1OverSqrt2);
    meshParmFile.write('R2OverSqrt2    %.12f;\n' % R2OverSqrt2);
    meshParmFile.write('R3OverSqrt2    %.12f;\n' % R3OverSqrt2);
    meshParmFile.write('z1    %.12f;\n' % z1);
    meshParmFile.write('z2    %.12f;\n' % z2);
    meshParmFile.write('nx    %d;\n' % nx);
    meshParmFile.write('ny    %d;\n' % ny);
    meshParmFile.write('nz    %d;\n' % nz);
    meshParmFile.write('z1PlusEps    %.12f;\n' % z1PlusEps);
    meshParmFile.write('z1MinusEps    %.12f;\n' % z1MinusEps);
    meshParmFile.write('z2PlusEps    %.12f;\n' % z2PlusEps);
    meshParmFile.write('z2MinusEps    %.12f;\n' % z2MinusEps);
    meshParmFile.write('R1PlusEps    %.12f;\n' % R1PlusEps);
    meshParmFile.write('mR1PlusEps    %.12f;\n' % mR1PlusEps);
    meshParmFile.write('grading    %.12f;\n' % grading);
