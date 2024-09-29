/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024
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
    Foam::floaterMotionRestraints::moorDynR2

Description
    floaterMotionRestraints model.  Mooring dynamics resolved.

    Calling MoorDyn: a lumped-mass based mooring dynamics model by
    M. Hall and A. Goupee, “Validation of a lumped-mass mooring line model with
    DeepCwind semisubmersible model test data,” Ocean Engineering, vol. 104,
    pp. 590–603, Aug. 2015.

    Chen, H., & Hall, M. (2022). CFD simulation of floating body motion with
    mooring dynamics: Coupling MoorDyn with OpenFOAM. Applied Ocean Research,
    124, 103210. https://doi.org/10.1016/j.apor.2022.103210

    Roenby, J., Aliyar, S., Bredmose, H. (2024). A robust algorithm for
    computational floating body dynamics. Royal Society Open Science,
    11:231453, https://doi.org/10.1098/rsos.231453

    Almeida Medina, T., Cercos-Pita, J.L., Gosset, A., Díaz Casás, V. (2024).
    A methodology for CFD simulations of floating structures in waves coupled
    with mooring dynamics.

SourceFiles
    moorDynR2.C

\*---------------------------------------------------------------------------*/

#ifndef moorDynR2_H
#define moorDynR2_H

#include "Time.H"
#include "pointField.H"
#include "floaterMotionRestraint.H"
#include "IFstream.H"
#include "OFstream.H"

#include "MoorDyn2.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace floaterMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class moorDynR2 Declaration
\*---------------------------------------------------------------------------*/

class moorDynR2
:
    public floaterMotionRestraint
{
    // Private data
        
        //- MoorDyn input filename
        fileName inputFile_;
        
        //- Coupling mode: 'POINT' or 'BODY'
        word couplingMode_;

        //- Number of coupling Dofs
        // 'POINT': nCouplingDof_ = 3 x npoints
        // 'BODY': nCouplingDof_ = 6 (1 body only)
        label nCouplingDof_;

        //- Reference point of attachment to the solid body
        pointField refAttachmentPt_;

        //- Initialize MoorDyn at the first time step ()
        mutable bool initialized_;
        
        //-  Write mooring lines config to file
        bool writeVTK_;
        
        //-  Format of the visualization file
        bool legacyVTK_;

        //- Mooring vtk filename prefix
        word vtkPrefix_;

        //- Start time to write VTK, default to runTime.start
        scalar vtkStartTime_;

        //- Sequence number to be included in vtk filename
        mutable label vtkCounter_;

        //- Current time, to count iteration
        mutable scalar curTime_;

        //- Iteration count, write VTK for the last iteration only
        mutable label iteration_;
        
        //- Outer corrector, times mooring states updated each time step
        label outerCorrector_;
        
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
            //moordyn_backup_.data = (uint64_t*)malloc(data_size);
            moordyn_backup_.data = static_cast<uint64_t*>(malloc(data_size));
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

        //- PVD file name
        inline fileName name_pvd() const
        {
            return "Mooring/VTK/" + vtkPrefix_ + ".pvd";
        }

        //- Check if a PVD file exists
        inline bool has_pvd() const
        {
            IFstream is_file(name_pvd());
            return is_file.good();
        }

        //- Create an empty PVD file
        inline void make_pvd() const
        {
            OFstream os(name_pvd());
            os << "<?xml version=\"1.0\"?>" << nl
            << "<VTKFile type=\"Collection\" version=\"0.1\"" << nl
            << "         byte_order=\"LittleEndian\"" << nl
            << "         compressor=\"vtkZLibDataCompressor\">" << nl
            << "  <Collection>" << nl
            << "  </Collection>" << nl
            << "</VTKFile>" << nl;
        }

        //- Read the PVD up to the part where new data whall be appended
        inline std::vector<string> read_pvd(const Time& time) const
        {
            IFstream is(name_pvd());
            if (!is.good()) {
                FatalError << "Error updating the PVD file: "
                    << "'moorDynR2.pvd' cannot be read"
                    << exit(FatalError);
            }
            std::vector<string> lines;
            while (true) {
                string line;
                is.getLine(line);
                if (line.ends_with("</Collection>")) {
                    return lines;
                }
                // Let's try to read the time to discard "the future"
                if (line.find("<DataSet") != std::string::npos) {
                    string stripped = line;
                    stripped.removeStart("    <DataSet timestep=\"");
                    stripped = stripped.substr(0, line.find("\""));
                    scalar t = std::stod(stripped);
                    if (t >= time.value()) {
                        return lines;
                    }
                }
                lines.push_back(line);
            }
            FatalError << "Error updating the PVD file: "
                << "'moorDynR2.pvd' is ill-formatted"
                << exit(FatalError);
        }

public:

    //- Runtime type information
    TypeName("moorDynR2");


    // Constructors

        //- Construct from components
        moorDynR2
        (
            const word& name,
            const dictionary& rBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<floaterMotionRestraint> clone() const
        {
            return autoPtr<floaterMotionRestraint>
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
            const floaterMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& rBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;
        
        //- Output vtk files for mooring
        void writeVTK(const Time&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
