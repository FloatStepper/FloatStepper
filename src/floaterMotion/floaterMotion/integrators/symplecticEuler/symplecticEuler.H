/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef symplecticEuler_H
#define symplecticEuler_H

#include "floaterIntegrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class symplecticEuler Declaration
\*---------------------------------------------------------------------------*/

class symplecticEuler : public floaterIntegrator
{
public:

    //- Runtime type information
    TypeName("symplecticEuler");

    // Constructors

        //- Construct from name and dictionary
        symplecticEuler(const word& name, const dictionary& dict);

        //- Construct as copy
        symplecticEuler(const symplecticEuler&) = default;

    //- Destructor
    virtual ~symplecticEuler() = default;

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
