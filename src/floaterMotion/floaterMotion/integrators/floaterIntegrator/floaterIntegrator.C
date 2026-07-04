/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#include "floaterIntegrator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(floaterIntegrator, 0);
    defineRunTimeSelectionTable(floaterIntegrator, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::floaterIntegrator::floaterIntegrator
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    coeffs_(dict.optionalSubDict(name + "Coeffs"))
{}

// ************************************************************************* //
