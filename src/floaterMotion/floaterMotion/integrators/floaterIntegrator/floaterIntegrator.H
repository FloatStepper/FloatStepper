/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef floaterIntegrator_H
#define floaterIntegrator_H

#include "floaterMotionState.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class floaterIntegrator Declaration
\*---------------------------------------------------------------------------*/

class floaterIntegrator
{
protected:

    // Protected data

        //- Name of the integrator
        const word name_;

        //- Reference to coefficients dictionary
        const dictionary& coeffs_;

public:

    //- Runtime type information
    TypeName("floaterIntegrator");

    declareRunTimeSelectionTable
    (
        autoPtr,
        floaterIntegrator,
        dictionary,
        (const word& name, const dictionary& dict),
        (name, dict)
    );

    // Constructors

        //- Construct from name and dictionary
        floaterIntegrator(const word& name, const dictionary& dict);

        //- Construct as copy
        floaterIntegrator(const floaterIntegrator&) = default;

    //- Destructor
    virtual ~floaterIntegrator() = default;

    // Static methods

        //- Return a reference to the selected integrator
        static autoPtr<floaterIntegrator> New
        (
            const word& integratorType,
            const dictionary& dict
        );

    // Member Functions

        //- Integrate the motion state given accelerations and time step
        virtual void integrate
        (
            floaterMotionState& motionState,
            const scalarField& dvwdt,
            const scalar deltaT
        ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
