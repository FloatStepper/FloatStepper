/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef basic_H
#define basic_H

#include "floaterIntegrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class basic Declaration
\*---------------------------------------------------------------------------*/

class basic : public floaterIntegrator
{
public:

    //- Runtime type information
    TypeName("basic");

    // Constructors

        //- Construct from name and dictionary
        basic(const word& name, const dictionary& dict);

        //- Construct as copy
        basic(const basic&) = default;

    //- Destructor
    virtual ~basic() = default;

    // Member Functions

        //- Integrate the motion state given accelerations and time step
        virtual void integrate
        (
            floaterMotionState& motionState,
            const scalarField& dvwdt,
            const scalar deltaT
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
