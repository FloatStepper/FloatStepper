/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef exponentialMap_H
#define exponentialMap_H

#include "floaterIntegrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class exponentialMap Declaration
\*---------------------------------------------------------------------------*/

class exponentialMap : public floaterIntegrator
{
public:

    //- Runtime type information
    TypeName("exponentialMap");

    // Constructors

        //- Construct from name and dictionary
        exponentialMap(const word& name, const dictionary& dict);

        //- Construct as copy
        exponentialMap(const exponentialMap&) = default;

    //- Destructor
    virtual ~exponentialMap() = default;

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
