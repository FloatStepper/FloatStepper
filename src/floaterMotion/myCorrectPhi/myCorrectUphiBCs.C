/*--------------------------------------------------- -----------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*--------------------------------------------------- -----------------------*/

#include "myCorrectPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::myCorrectUphiBCs
(
    volVectorField& U,
    surfaceScalarField& phi
)
{
    const fvMesh& mesh = U.mesh();

    if (mesh.changing())
    {
        volVectorField::Boundary& Ubf = U.boundaryFieldRef();
        surfaceScalarField::Boundary& phibf =
            phi.boundaryFieldRef();

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                Ubf[patchi].initEvaluate();
            }
        }

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                Ubf[patchi].evaluate();

                phibf[patchi] = Ubf[patchi] & mesh.Sf().boundaryField()[patchi];
            }
        }
    }
}


void Foam::myCorrectUphiBCs
(
    const volScalarField& rho,
    volVectorField& U,
    surfaceScalarField& phi
)
{
    const fvMesh& mesh = U.mesh();

    if (mesh.changing())
    {
        volVectorField::Boundary& Ubf = U.boundaryFieldRef();
        surfaceScalarField::Boundary& phibf =
            phi.boundaryFieldRef();

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                Ubf[patchi].initEvaluate();
            }
        }

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                Ubf[patchi].evaluate();

                phibf[patchi] =
                    rho.boundaryField()[patchi]
                   *(
                        Ubf[patchi]
                      & mesh.Sf().boundaryField()[patchi]
                    );
            }
        }
    }
}


// ************************************************************************* //
