/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#include "myCorrectPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvScalarMatrix.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "myAdjustPhi.H"
#include "fvcMeshPhi.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RAUfType, class DivUType>
void Foam::myCorrectPhi
(
    volVectorField& U,
    surfaceScalarField& phi,
    const volScalarField& p,
    const RAUfType& rAUf,
    const DivUType& divU,
    pimpleControl& pimple
)
{
    const fvMesh& mesh = U.mesh();
    const Time& runTime = mesh.time();

    myCorrectUphiBCs(U, phi);

    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
        // The line above is replaced by the outcommented line below from
        // v2312, but this prevents compilation with v2206-v2306
//        fvPatchFieldBase::zeroGradientType()
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions(), Zero),
        pcorrTypes
    );

    if (pcorr.needReference())
    {
        fvc::makeRelative(phi, U);
        myAdjustPhi(phi, U, pcorr);
        fvc::makeAbsolute(phi, U);
    }

    mesh.setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi) - divU
        );

        pcorrEqn.setReference(0, 0);

        pcorrEqn.solve
        (
            mesh.solver(pcorr.select(pimple.finalNonOrthogonalIter()))
            // The line above is replaced by the outcommented line below from
            // v2312, but this prevents compilation with v2206-v2306
//            pcorr.select(pimple.finalNonOrthogonalIter())
);

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }
}


template<class RAUfType, class DivRhoUType>
void Foam::myCorrectPhi
(
    volVectorField& U,
    surfaceScalarField& phi,
    const volScalarField& p,
    const volScalarField& rho,
    const volScalarField& psi,
    const RAUfType& rAUf,
    const DivRhoUType& divRhoU,
    pimpleControl& pimple
)
{
    const fvMesh& mesh = U.mesh();
    const Time& runTime = mesh.time();

    myCorrectUphiBCs(rho, U, phi);

    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
        // The line above is replaced by the outcommented line below from
        // v2312, but this prevents compilation with v2206-v2306
//        fvPatchFieldBase::zeroGradientType()
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions(), Zero),
        pcorrTypes
    );

    mesh.setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divRhoU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::ddt(psi, pcorr)
          + fvc::div(phi)
          - fvm::laplacian(rAUf, pcorr)
         ==
            divRhoU
        );

        pcorrEqn.solve
        (
            mesh.solver(pcorr.select(pimple.finalNonOrthogonalIter()))
            // The line above is replaced by the outcommented line below from
            // v2312, but this prevents compilation with v2206-v2306
//            pcorr.select(pimple.finalNonOrthogonalIter())
);

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pcorrEqn.flux();
        }
    }
}


// ************************************************************************* //
