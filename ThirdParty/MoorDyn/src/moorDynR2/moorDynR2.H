/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef moorDynR2_H
#define moorDynR2_H

#include "rigidBodyMotionRestraint.H"
#include "MoorDyn2.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace rigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class moorDynR2 Declaration
\*---------------------------------------------------------------------------*/

class moorDynR2
:
    public rigidBodyMotionRestraint
{
    // Private data

	//- Initialize MoorDyn at the first time step ()
        mutable bool initialized_;


public:

    //- Runtime type information
    TypeName("moorDynR2");


    // Constructors

        //- Construct from components
        moorDynR2
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<rigidBodyMotionRestraint> clone() const
        {
            return autoPtr<rigidBodyMotionRestraint>
            (
                new moorDynR2(*this)
            );
        }


    //- Destructor
    virtual ~moorDynR2();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const rigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;

private:

    // The MoorDyn system instance
    MoorDyn moordyn_;
    // The coupled body to be controlled by this MoorDyn instance
    MoorDynBody moordyn_body_;
    // A backup of the Mooring system to can restore in between implicit time
    // steps
    typedef struct _moordyn_backup {
        double t;
        uint64_t *data;
    } moordyn_backup;
    mutable moordyn_backup moordyn_backup_;

    //- Create a backup of the mooring system
    inline void save_mooring(const double t) const
    {
        int err;
        size_t data_size;
        err = MoorDyn_Serialize(moordyn_, &data_size, NULL);
        if (err != MOORDYN_SUCCESS) {
            FatalError << "MoorDyn system size cannot be retrieved"
                       << exit(FatalError);
        }
        if (moordyn_backup_.data)
            free(moordyn_backup_.data);
        moordyn_backup_.data = (uint64_t*)malloc(data_size);
        if (!moordyn_backup_.data) {
            FatalError << "Failure allocating " << data_size << " bytes"
                       << exit(FatalError);
        }
        err = MoorDyn_Serialize(moordyn_, NULL, moordyn_backup_.data);
        if (err != MOORDYN_SUCCESS) {
            FatalError << "MoorDyn system cannot be saved"
                       << exit(FatalError);
        }
        moordyn_backup_.t = t;
    }

    //- Restore the backup of the mooring system
    inline void load_mooring() const
    {
        if (!moordyn_backup_.data) {
            FatalError << "No available backup"
                       << exit(FatalError);
        }
        const int err = MoorDyn_Deserialize(moordyn_, moordyn_backup_.data);
        if (err != MOORDYN_SUCCESS) {
            FatalError << "MoorDyn system cannot be restored"
                       << exit(FatalError);
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
