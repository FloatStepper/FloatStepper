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

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
References

    Chen, H., & Hall, M. (2022). CFD simulation of floating body motion with
    mooring dynamics: Coupling MoorDyn with OpenFOAM. Applied Ocean Research,
    124, 103210. https://doi.org/10.1016/j.apor.2022.103210

\*---------------------------------------------------------------------------*/

#include "moorDynR2.H"
#include "addToRunTimeSelectionTable.H"
#include "floaterMotion.H"
//#include "sixDoFRigidBodyMotion.H"
//#include "Time.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "error.H"
#include "quaternion.H"
#include "foamVtkSeriesWriter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace floaterMotionRestraints
{
    defineTypeNameAndDebug(moorDynR2, 0);

    addToRunTimeSelectionTable
    (
        floaterMotionRestraint,
        moorDynR2,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::moorDynR2::moorDynR2
(
    const word& name,
    const dictionary& rBMRDict
)
:
    floaterMotionRestraint(name, rBMRDict)
{
    read(rBMRDict);
    
    initialized_ = false;
    vtkCounter_ = 0;
    curTime_ = -1;
    iteration_ = 0;

    Info<< "Create moorDynR2 using MoorDyn v2."
        << "  Coupling mode: " << couplingMode_ 
        << endl;
    
    moordyn_backup_.t = 0.0;
    moordyn_backup_.data = nullptr;

    // Create the MoorDyn system
    if (Pstream::master())
    {
        int moordyn_err = MOORDYN_SUCCESS;
        //moordyn_ = MoorDyn_Create("Mooring/lines_v2.txt");
        moordyn_ = MoorDyn_Create(inputFile_.c_str());
        if (!moordyn_) {
            FatalError << "MoorDyn v2 cannot be created!" << exit(FatalError);
        }
        // In this release just a single floating body is accepted. In the
        // future more bodies can be added, providing info about the body to
        // be modelled
        unsigned int n;
        moordyn_err = MoorDyn_NCoupledDOF(moordyn_, &n);
        if (moordyn_err != MOORDYN_SUCCESS) {
            FatalErrorInFunction
                << "NCoupledDOF error on the MoorDyn definition file" 
                << exit(FatalError);
        }
        nCouplingDof_ = int(n);

        if (couplingMode_ == word("POINT"))
        {
            Info<< "\tCoupling mode: " << couplingMode_ << ", expecting nCouplingDof="
                << 3 * refAttachmentPt_.size() << endl;
            if (nCouplingDof_ != 3 * refAttachmentPt_.size())
            {
                FatalErrorInFunction
                    << " refAttachmentPt size: " << refAttachmentPt_.size() 
                    << " not equal to # coupled points " << int(nCouplingDof_/3)
                    << " defined in MoorDyn input!"
                    << exit(FatalError);
            }
        }
        else if (couplingMode_ == word("BODY"))
        {
            Info<< "Coupling mode: " << couplingMode_ << ", expecting nCouplingDof = 6"
                << endl;
            if (nCouplingDof_ != 6 )
            {
                FatalErrorInFunction
                    << " Support one coupled only in this sixDoFMotion restraint"
                    << exit(FatalError);
            }

            moordyn_err = MoorDyn_GetNumberBodies(moordyn_, &n);
            if ((moordyn_err != MOORDYN_SUCCESS) || n != 1) {
                FatalError << "Only one coupled body was expected "
                        << "on the MoorDyn definition file" << exit(FatalError);
            }
            moordyn_body_ = NULL;
            for (unsigned int i = 0; i < n; i++) {
                MoorDynBody body = MoorDyn_GetBody(moordyn_, i + 1);
                if (!body) {
                    FatalError << "Failure getting the MoorDyn body " << i + 1
                            << exit(FatalError);
                }
                int t;
                moordyn_err = MoorDyn_GetBodyType(body, &t);
                if (moordyn_err != MOORDYN_SUCCESS) {
                    FatalError << "Failure geeting the body " << i + 1
                            << " type" << exit(FatalError);
                }
                if (t == -1) {
                    // Coupled body, see Body.hpp:143
                    moordyn_body_ = body;
                    break;
                }
            }
            if (!moordyn_body_) {
                FatalError << "No coupled body could be found" << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Two coupling modes supported: 'POINT' or 'BODY' "
                << exit(FatalError);
        }   

        if (writeVTK_)
        {
#ifndef MOORDYN2_HAVE_VTK
            if (!legacyVTK_) {
                FatalErrorInFunction
                    << "vtkLegacyFormat=false option has been set. "
                    << "However, MoorDyn v2 has been compiled without VTK "
                    << "support. Please set vtkLegacyFormat=true or install "
                    << "VTK and recompile foamMooring"
                    << exit(FatalError);
            }
#endif
        }

    }

    if (writeVTK_)
    {
        mkDir("Mooring/VTK");
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::moorDynR2::~moorDynR2()
{
    if (initialized_)
    {
        // Close MoorDyn call
        MoorDyn_Close(moordyn_);
        if (moordyn_backup_.data)
            free(moordyn_backup_.data);

        vtk::seriesWriter writer;

        writer.scan("Mooring/VTK/" + vtkPrefix_ + ".vtk");
        
        Info<< "Writing mooring vtk series file" << nl;

        fileName vtkSeries("Mooring/VTK/" + vtkPrefix_ + ".vtk");
        writer.write(vtkSeries);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::floaterMotionRestraints::moorDynR2::restrain
(
    const floaterMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    if (Pstream::master())
{
    int moordyn_err = MOORDYN_SUCCESS;
    
    const Time& time = motion.time();

    scalar deltaT = time.deltaTValue();
    scalar t = time.value();
    scalar tprev = t - deltaT;
    
    if (t == curTime_)
    {
        iteration_++;
    }
    else if (t > curTime_)
    {
        iteration_ = 1;
        curTime_ = t;
    }

    pointField fairPos = vectorField(int(nCouplingDof_/3), vector::zero);
    vectorField fairVel = vectorField(int(nCouplingDof_/3), vector::zero);
    vectorField fairForce =vectorField(int(nCouplingDof_/3), vector::zero);

    // If coupling mode is 'BODY', X and fairPos are in fact body's 6DoF motion.
    // Flines is mooring forces and moments on the body (may need to reverse the sign).
    double* X = &fairPos[0][0];
    double* XD = &fairVel[0][0];
    double* Flines = &fairForce[0][0];

    if (couplingMode_ == word("BODY"))
    {
        point CoM = motion.centreOfMass();
        vector rotationAngle
        (
            quaternion(motion.orientation()).eulerAngles(quaternion::XYZ)
        );

        vector v = motion.v();
        vector omega = motion.omega();

        for (int ii=0; ii<3; ii++)
        {
            X[ii] = CoM[ii];
            X[ii+3] = rotationAngle[ii];
            XD[ii] = v[ii];
            XD[ii+3] = omega[ii];
        }
    }
    else
    {
        // Calculate fairlead position
        for(int pt=0; pt<refAttachmentPt_.size(); pt++)
        {
            fairPos[pt] = motion.transform(refAttachmentPt_[pt]);
            fairVel[pt] = motion.velocity(refAttachmentPt_[pt]);
        }
    }

    Info<< "\n\ttprev = " << tprev << ", X[6]/fairPosition: " << fairPos << endl;

    if (!initialized_)
    {
        moordyn_err = MoorDyn_Init(moordyn_, X, XD);
        if (moordyn_err != MOORDYN_SUCCESS) {
            FatalError << "MoorDyn could not be initialized"
                       << exit(FatalError);
        }
        Info<< "MoorDyn module initialized!" << endl;
        initialized_ = true;
        save_mooring(tprev);
        if (writeVTK_) {
            autoPtr<Time> dummy_t = Time::New();
            dummy_t->setTime(0.0, 0);
            writeVTK(dummy_t.ref());
        }

        curTime_ = t;
        iteration_ = 1;
    } else if (tprev - moordyn_backup_.t >= 1.e-3 * deltaT) {
        // We have successfully advanced forward in time
        save_mooring(tprev);
        Info<< "MoorDyn module saved at t = " << tprev << " s" << endl;
    } else {
        // We are repeating the same time step because the implicit scheme
        load_mooring();
        Info<< "MoorDyn module restored to t = " << moordyn_backup_.t << " s" << endl;
    }
    
    // Step MoorDyn to get mooring forces on body
    //moordyn_err = MoorDyn_Step(moordyn_, X, XD, Flines, &tprev, &deltaT);
    moordyn_err = MoorDyn_Step(moordyn_, &fairPos[0][0], &fairVel[0][0], &fairForce[0][0], &tprev, &deltaT);
    if (moordyn_err != MOORDYN_SUCCESS) {
        FatalErrorInFunction
            << "Error computing MoorDyn step " << tprev
            << "s -> " << tprev + deltaT << "s" << exit(FatalError);
    }

    if (couplingMode_ == word("BODY"))
    {
        for(int i=0;i<3;i++)
        {
            restraintForce[i] = Flines[i];
            restraintMoment[i] = Flines[i+3];
        }
    }
    else
    {
        // Sum forces and calculate moments
        point CoR = motion.centreOfRotation();

        for(int pt=0; pt<refAttachmentPt_.size(); pt++)
        {
            restraintForce += fairForce[pt];

            restraintMoment += (fairPos[pt] - CoR) ^ fairForce[pt];

            Info<< " attachPt[" << pt << "] " // << fairPos[pt]
                << " force " << fairForce[pt]
                //<< " moment " << moment
                << endl; 
        }
    }

    // Since moment is already calculated by LinesCalc, set to
    // centreOfRotation to be sure of no spurious moment
    restraintPosition = motion.centreOfRotation();

    //if (motion.report())
    {
        Info << t << ": force " << restraintForce
             << ", moment " << restraintMoment
             << endl;
    }
    
    if (writeVTK_ && iteration_ == outerCorrector_)
    {
        if (t >= vtkStartTime_ && time.outputTime())
        {
            Info<< "Write mooring VTK ..." << endl;
            writeVTK(time);
        }
    }
}
    //Distribute results to other nodes:
    Pstream::broadcast(restraintPosition);
    Pstream::broadcast(restraintForce);
    Pstream::broadcast(restraintMoment);
}


bool Foam::floaterMotionRestraints::moorDynR2::read
(
    const dictionary& rBMRDict
)
{
    floaterMotionRestraint::read(rBMRDict);
    
    rBMRCoeffs_.readEntry("inputFile", inputFile_);
    rBMRCoeffs_.readEntry("couplingMode", couplingMode_);
    if (couplingMode_ == word("POINT"))
    {
        rBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
    }

    writeVTK_ = rBMRCoeffs_.getOrDefault<Switch>("writeMooringVTK", false);
    if (writeVTK_)
    {
#ifdef MOORDYN2_HAVE_VTK
        const bool default_legacy_vtk = false;
#else
        const bool default_legacy_vtk = true;
#endif
        vtkPrefix_ = rBMRCoeffs_.getOrDefault<word>("vtkPrefix", "mdV2");
        vtkStartTime_ = rBMRCoeffs_.getOrDefault<scalar>("vtkStartTime", 0);
        legacyVTK_ = rBMRCoeffs_.getOrDefault<Switch>("vtkLegacyFormat",
                                                          default_legacy_vtk);
        outerCorrector_ = rBMRCoeffs_.getOrDefault<scalar>("outerCorrector", 1);
    }
    
    return true;
}


void Foam::floaterMotionRestraints::moorDynR2::write
(
    Ostream& os
) const
{
    os.writeEntry("inputFile", inputFile_);
    os.writeEntry("couplingMode", couplingMode_);
    os.writeEntry("nCouplingDof", nCouplingDof_);
    if (couplingMode_ == word("POINT"))
        os.writeEntry("refAttachmentPt", refAttachmentPt_);
    
    os.writeEntry("writeMooringVTK", writeVTK_);
    if (writeVTK_)
    {
        os.writeEntry("vtkPrefix", vtkPrefix_);
        os.writeEntry("vtkStartTime", vtkStartTime_);
        os.writeEntry("vtkLegacyFormat", legacyVTK_);
        os.writeEntry("outerCorrector", outerCorrector_);
    }
}

void Foam::floaterMotionRestraints::moorDynR2::writeVTK(const Time& time) const
{
    fileName name(
        vtkPrefix_ + "_" + Foam::name(++vtkCounter_));
    if (!legacyVTK_) {
        name += ".vtm";
        int moordyn_err = MoorDyn_SaveVTK(moordyn_,
                                          ("Mooring/VTK/" + name).c_str());
        if (moordyn_err != MOORDYN_SUCCESS) {
            FatalError << "Error saving the VTK file \""
                       << name << " for time " << time.timeName() << " s"
                       << exit(FatalError);
        }
        // Update the vtp file
        if (!has_pvd()) {
            make_pvd();
        }
        auto lines = read_pvd(time);
        OFstream os(name_pvd());
        for (auto line : lines) {
            os << line.c_str() << nl;
        }
        os << "    <DataSet timestep=\"" << time.timeName()
           << "\" group=\"\" part=\"0\" "
           << "file=\"" << name.c_str() << "\"/>" << nl
           << "  </Collection>" << nl
           << "</VTKFile>" << nl;

        return;
    }

    name += ".vtk";
    OFstream mps("Mooring/VTK/" + name);
    mps.precision(4);

    unsigned int nLines, nSeg;
    MoorDyn_GetNumberLines(moordyn_, &nLines);

    labelList nodesPerLine(nLines, -1);
    for(int i=0; i<int(nLines); i++)
    {
        MoorDynLine line = MoorDyn_GetLine(moordyn_, i+1);
        MoorDyn_GetLineN(line, &nSeg);
        nodesPerLine[i] = nSeg+1;
    }

    double coord[max(nodesPerLine)][3];

    // Writing header
    mps << "# vtk DataFile Version 3.0" << nl
        << "MoorDyn v2 vtk output time=" << time.timeName()
        << nl << "ASCII" << nl << "DATASET POLYDATA" << endl;
 
    // Writing points
    mps << "\nPOINTS " << sum(nodesPerLine) << " float" << endl;

    for(int i=0; i<int(nLines); i++)
    {   
        //map_.getNodeCoordinates(i, nodesPerLine_[i], &coord[0][0]);
        MoorDynLine line = MoorDyn_GetLine(moordyn_, i+1);
        
        for(int p=0; p<nodesPerLine[i]; p++)
        {
            MoorDyn_GetLineNodePos(line, p, &coord[p][0]);
            
            mps << coord[p][0] << " " << coord[p][1] << " " << coord[p][2] << endl;
        }
    }
    
    // Writing lines
    mps << "\nLINES " << nLines << " " << sum(nodesPerLine+1) << endl;

    label start_node(0);
    for(int i=0; i<int(nLines); i++)
    {       
        mps << nodesPerLine[i];
        
        for(int j=0; j<nodesPerLine[i]; j++)
        {
            mps << " " << start_node+j;
        }
        mps << endl;
        
        start_node += nodesPerLine[i];
    }
}

// ************************************************************************* //
