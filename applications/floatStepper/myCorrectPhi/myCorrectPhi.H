/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef myCorrectPhi_H
#define myCorrectPhi_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class pimpleControl;

    //- If the mesh is moving correct the velocity BCs on the moving walls to
    //  ensure the corrected fluxes and velocity are consistent
    void myCorrectUphiBCs
    (
        volVectorField& U,
        surfaceScalarField& phi
    );

    //- If the mesh is moving correct the velocity BCs on the moving walls to
    //  ensure the corrected fluxes and velocity are consistent
    void myCorrectUphiBCs
    (
        const volScalarField& rho,
        volVectorField& U,
        surfaceScalarField& phi
    );

    template<class RAUfType, class DivUType>
    void myCorrectPhi
    (
        volVectorField& U,
        surfaceScalarField& phi,
        const volScalarField& p,
        const RAUfType& rAUf,
        const DivUType& divU,
        pimpleControl& pimple
    );

    template<class RAUfType, class DivRhoUType>
    void myCorrectPhi
    (
        volVectorField& U,
        surfaceScalarField& phi,
        const volScalarField& p,
        const volScalarField& rho,
        const volScalarField& psi,
        const RAUfType& rAUf,
        const DivRhoUType& divRhoU,
        pimpleControl& pimple
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myCorrectPhi.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
