/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef morphingBoxRegion_H
#define morphingBoxRegion_H

#include "coordinateSystem.H"
#include "dictionary.H"
#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Struct boxBounds Declaration
\*---------------------------------------------------------------------------*/

struct boxBounds
{
    scalar xMin, xMax;
    scalar yMin, yMax;
    scalar zMin, zMax;
};


/*---------------------------------------------------------------------------*\
                    Class morphingBoxRegion Declaration
\*---------------------------------------------------------------------------*/

class morphingBoxRegion
{
    // Private data

        //- Coordinate system for this morphing region
        coordinateSystem cs_;

        //- Whether bounds are relative to "body" or "lab" frame
        word relativeTo_;

        //- Inner box bounds
        boxBounds inner_;

        //- Outer box bounds
        boxBounds outer_;


    // Private member functions

        //- Read a boxBounds from a sub-dictionary
        static boxBounds readBounds
        (
            const dictionary& dict,
            const scalar defaultValue
        );


public:

    // Constructors

        //- Construct from dictionary, with optional fallback dictionary
        //  for entries not present in dict (implements $xMorphing inheritance)
        morphingBoxRegion
        (
            const dictionary& dict,
            const dictionary* fallbackPtr = nullptr
        );


    // Member Functions

        //- Update a scalar blending field in place.
        //  s is multiplied by the blending weight derived from all six
        //  faces of the morphing box. Points are assumed to already be
        //  expressed in the frame consistent with relativeTo.
        void updateBlending
        (
            const pointField& points,
            scalarField& s
        ) const;

        //- Shift boxes e.g. if coordinates are relative to body position
        void shiftBoxes(const vector& cor);

    // Access

        //- Coordinate system
        const coordinateSystem& cs() const { return cs_; }

        //- Reference frame: "body" or "lab"
        const word& relativeTo() const { return relativeTo_; }

        //- Inner box bounds
        const boxBounds& inner() const { return inner_; }

        //- Outer box bounds
        const boxBounds& outer() const { return outer_; }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
