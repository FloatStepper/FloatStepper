/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Solver for the ODEs governing the motion of an ellipse shaped rigid body
    moving in an infinite, ideal 2D fluid. These equations are sometimes called 
    the Kirchhoff-Kelvin equations. The degrees of freedom are the two 
    components of the body geometric centre and its orientation. The body may
    be subject to a constant gravity and buoyancy force.

    The input parameters are provided by a bodyDict an exampe of which is 
    provided in the solver directory.

    For the meaning of the various parameters please see the comments in the 
    code below.

Author
    Johan Roenby, STROMNING, 2021

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "diagTensor.H"
#include "mathematicalConstants.H"
#include "Fstream.H"

using namespace Foam;
using namespace Foam::constant::mathematical;
using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class KirchhoffKelvin2D
:
    public ODESystem
{
    //- Added mass diagonal coefficients
    scalar A11_;
    scalar A22_;
    scalar A33_;

    //- Body mass
    scalar Mb_;

    //- Body moment of inertia relative to centre
    scalar Icm_;

    //- Effective mass diagonal coefficients
    scalar M11_;
    scalar M22_;
    scalar M33_;

    //- Net force along y-axis (buoyancy)
    scalar Fy_;

public:

    KirchhoffKelvin2D()
    {}

    KirchhoffKelvin2D
    (
        dictionary& dict
    )
    {
        // Body radius
        scalar R(dict.get<scalar>("R"));

        // Eccentricity (in complex notation z = R(zeta + a^2/zeta)
        scalar a(dict.get<scalar>("a"));

        // Fluid density
        scalar rhof(dict.get<scalar>("rhof"));

        // Adding added mass coefficients
        A11_ = rhof*pi*Foam::sqr(R*(1.0 - a*a));
        A22_ = rhof*pi*Foam::sqr(R*(1.0 + a*a));
        A33_ = 2.0*rhof*pi*Foam::pow4(R*a);

        // Option to specify added mass coefficients
        A11_ = dict.getOrDefault<scalar>("A11", A11_);
        A22_ = dict.getOrDefault<scalar>("A22", A22_);
        A33_ = dict.getOrDefault<scalar>("A33", A33_);

        Info << "Added mass: A11 = " << A11_ << ", A22 = " << A22_
            << ", A33 = " << A33_ << endl;

        // Gravity component along y-axis (negative = downwards)
        scalar gy(dict.getOrDefault<scalar>("gy", 0));
        // Body mass density (assumed uniform)
        scalar rhob(dict.get<scalar>("rhob"));
        // Body area
        scalar Ab = pi*R*R*(1.0 - pow4(a));
        // Body mass
        Mb_ = rhob*Ab;
        // Body moment of inertia
        Icm_ = 0.5*rhob*pi*pow4(R)*(1 - pow(a,8));

        // Option to specify body mass and moment of inertia
        Mb_ = dict.getOrDefault<scalar>("Mb", Mb_);
        Icm_ = dict.getOrDefault<scalar>("Icm", Icm_);

        Info << "Body mass: Mb = " << Mb_ << endl;
        Info << "Body moment of inertia: I = " << Icm_ << endl;

        // Adding mass and inertia to efffective mass
        M11_ = A11_ + Mb_;
        M22_ = A22_ + Mb_;
        M33_ = A33_ + Icm_;

        // Possibility to specify effective mass - overwriting above values
        M11_ = dict.getOrDefault<scalar>("M11", M11_);
        M22_ = dict.getOrDefault<scalar>("M22", M22_);
        M33_ = dict.getOrDefault<scalar>("M33", M33_);

        //Buoyancy force
        Fy_ = Ab*(rhob - rhof)*gy;

        Info << "Effective mass: M11 = " << M11_ << ", M22 = " << M22_
            << ", M33 = " << M33_ << ", Fy = " << Fy_ << endl;

        scalarField Y(nEqns());
        Y[0] = dict.get<scalar>("x0");
        Y[1] = dict.get<scalar>("y0");
        Y[2] = pi*dict.get<scalar>("th0ByPi");
        Y[3] = dict.get<scalar>("vx0");
        Y[4] = dict.get<scalar>("vy0");
        Y[5] = dict.get<scalar>("omega0");

        scalar x0 = Y[0];
        scalar y0 = Y[1];
        scalar th0 = Y[2];
        scalar vx0 = Y[3];
        scalar vy0 = Y[4];
        scalar omega0 = Y[5];

        // Calculating initial accelerations (needed for CFD initialisation)
        scalarField dYdt0(nEqns(), 0);
        derivatives(0, Y, dYdt0);
        scalar ax0 = dYdt0[3];
        scalar ay0 = dYdt0[4];
        scalar alpha0 = dYdt0[5];

        // calculating initial body angular momentum and torque (needed by
        // sixDoFRigidBodyMotionState)
        scalar L = Icm_*omega0;
        scalar costh = Foam::cos(th0);
        scalar sinth = Foam::sin(th0);
        scalar U1 = vx0*costh + vy0*sinth;
        scalar U2 = -vx0*sinth + vy0*costh;
        scalar tau0 = (A11_ - A22_)*U1*U2 - A33_*alpha0;

        // Writing file for initialising CFD with rigidBodyMoion
        fileName rbmFileName("rigidBodyMotionStateInput");
        dictionary rbmDict(rbmFileName);
        rbmDict.add("centreOfRotation", vector(x0, y0, 0), true);
        rbmDict.add("orentation", tensor(costh, -sinth, 0, sinth, costh, 0, 0, 0, 1), true);
        rbmDict.add("velocity", vector(vx0, vy0, 0), true);
        rbmDict.add("acceleration", vector(ax0, ay0, 0), true);
        rbmDict.add("omega", vector(0, 0, omega0), true);
        rbmDict.add("domegadt", vector(0, 0, alpha0), true);
        rbmDict.write(OFstream(rbmFileName)(), false);

        // Writing file for initialising CFD with sixDoFRigidBodyMoion
        fileName sixDoFFileName("sixDoFRigidBodyMotionStateInput");
        dictionary sixDoFDict(sixDoFFileName);
        sixDoFDict.add("centreOfRotation", vector(x0, y0, 0), true);
        sixDoFDict.add("orentation", tensor(costh, -sinth, 0, sinth, costh, 0, 0, 0, 1), true);
        sixDoFDict.add("velocity", vector(vx0, vy0, 0), true);
        sixDoFDict.add("acceleration", vector(ax0, ay0, 0), true);
        sixDoFDict.add("angularMomentum", vector(0, 0, L), true);
        sixDoFDict.add("torque", vector(0, 0, tau0), true);
        sixDoFDict.write(OFstream(sixDoFFileName)(), false);

    }

    label nEqns() const
    {
        return 6;
    }

    void derivatives
    (
        const scalar t,
        const scalarField& Y,
        scalarField& dYdt
    ) const
    {
        // Equations from matlab script
        // Q = [cos(th) -sin(th); sin(th) cos(th)];
        // Q.T() = [cos(th) sin(th); -sin(th) cos(th)];
        // V = [v1; v2];
        // U = Q'*V;
        // G = Q'*F;
        // U = Q'V = [cos(th) sin(th); -sin(th) cos(th)] V;
        //dUdt1 = (A(2)*U(2)*omega + G(1))/A(1);
        //dUdt2 = (-A(1)*U(1)*omega + G(2))/A(2);
        //dUdt = [dUdt1; dUdt2];
        //dVdt = Q*dUdt - omega*[v2; -v1];

        //domegadt = ((A(1)-A(2))*U(1)*U(2))/A(3);

        // scalar x = Y[0];
        // scalar y = Y[1];
        scalar th = Y[2];
        scalar v1 = Y[3];
        scalar v2 = Y[4];
        scalar omega = Y[5];

        scalar costh = Foam::cos(th);
        scalar sinth = Foam::sin(th);

        // Velocity vector components in body-fixed coordinates
        scalar U1 = v1*costh + v2*sinth;
        scalar U2 = -v1*sinth + v2*costh;

        // Gravity components in body-fixed coordinates
        scalar Fx = 0.0;
        scalar G1 = Fx*costh + Fy_*sinth;
        scalar G2 = -Fx*sinth + Fy_*costh;

        // Linear acceleration components
        scalar dUdt1 = (M22_*U2*omega + G1)/M11_;
        scalar dUdt2 = (-M11_*U1*omega + G2)/M22_;
        scalar dv1dt = costh*dUdt1 - sinth*dUdt2 - omega*v2;
        scalar dv2dt = sinth*dUdt1 + costh*dUdt2 + omega*v1;

        // Angular acceleration
        scalar domegadt = ((M11_ - M22_)*U1*U2)/M33_;

        dYdt[0] = v1;
        dYdt[1] = v2;
        dYdt[2] = omega;
        dYdt[3] = dv1dt;
        dYdt[4] = dv2dt;
        dYdt[5] = domegadt;

    }

    void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {
        dfdx = 0.0;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addArgument("ODESolver");
    argList args(argc, argv);
    dictionary dict;
    dict.add("solver", args[1]);

    dictionary bodyDict(IFstream("bodyDict")());
    KirchhoffKelvin2D ode(bodyDict);

    // Create the selected ODE system solver
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict);

    odeSolver->relTol() = 1e-12;
    odeSolver->absTol() = 1e-12;
    scalar dtEst = 1e-8;

    // Initialise the ODE system fields
    scalarField Y(ode.nEqns());
    Y[0] = bodyDict.get<scalar>("x0");
    Y[1] = bodyDict.get<scalar>("y0");
    Y[2] = pi*bodyDict.get<scalar>("th0ByPi");
    Y[3] = bodyDict.get<scalar>("vx0");
    Y[4] = bodyDict.get<scalar>("vy0");
    Y[5] = bodyDict.get<scalar>("omega0");

    scalar tStart = bodyDict.getOrDefault<scalar>("tstart", 0);

    Info << "Solution: t = " << tStart << ", Y = " << Y << endl;

    scalar tEnd = bodyDict.get<scalar>("tend");
    scalar nPlotSteps = bodyDict.get<scalar>("nPlotSteps");

    OFstream outFile(bodyDict.get<word>("outFile"));

    outFile << "# Time x y th vx vy omega" << endl;
    scalar ti = tStart;
    for (label i=0; i<nPlotSteps; i++)
    {
        scalar tip1 = tStart + (i+1)*(tEnd-tStart)/nPlotSteps;
        odeSolver->solve(ti, tip1, Y, dtEst);
        Info << "Solution: t = " << tip1 << ", Y = " << Y << endl;
        outFile << tip1 << " " 
            << Y[0] << " "
            << Y[1] << " "
            << Y[2] << " "
            << Y[3] << " "
            << Y[4] << " "
            << Y[5] << " "
            << endl;
        ti = tip1;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
