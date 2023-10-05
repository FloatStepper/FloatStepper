/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::rigidBodyMotionRestraints::moorDynR2

Description
    rigidBodyMotionRestraints model.  Mooring dynamics resolved.

    Calling MoorDyn: a lumped-mass based mooring dynamics model by
    M. Hall and A. Goupee, “Validation of a lumped-mass mooring line model with
    DeepCwind semisubmersible model test data,” Ocean Engineering, vol. 104,
    pp. 590–603, Aug. 2015.

    Chen, H., & Hall, M. (2022). CFD simulation of floating body motion with mooring dynamics: Coupling 
    MoorDyn with OpenFOAM. Applied Ocean Research, 124, 103210.
    https://doi.org/10.1016/j.apor.2022.103210

SourceFiles
    moorDynR2.C

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
