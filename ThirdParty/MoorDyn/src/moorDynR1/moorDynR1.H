/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef moorDynR1_H
#define moorDynR1_H

#include "floaterMotionRestraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace floaterMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class moorDynR1 Declaration
\*---------------------------------------------------------------------------*/

class moorDynR1
:
    public floaterMotionRestraint
{
    // Private data
    vector initialCoM;

    vector DecayTranslation_;

    vector DecayRotation_;

    //- Initialize MoorDyn at the first time step ()
    mutable bool initialized_;


public:

    //- Runtime type information
    TypeName("moorDynR1");


    // Constructors

        //- Construct from components
        moorDynR1
        (
            const word& name,
            const dictionary& rBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<floaterMotionRestraint> clone() const
        {
            return autoPtr<floaterMotionRestraint>
            (
                new moorDynR1(*this)
            );
        }


    //- Destructor
    virtual ~moorDynR1();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const floaterMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& rBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
