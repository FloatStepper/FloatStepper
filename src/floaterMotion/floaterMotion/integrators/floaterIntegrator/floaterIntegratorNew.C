/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#include "floaterIntegrator.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::floaterIntegrator>
Foam::floaterIntegrator::New
(
    const word& integratorType,
    const dictionary& dict
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(integratorType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "floaterIntegrator",
            integratorType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<floaterIntegrator>(cstrIter()(integratorType, dict));
}

// ************************************************************************* //
