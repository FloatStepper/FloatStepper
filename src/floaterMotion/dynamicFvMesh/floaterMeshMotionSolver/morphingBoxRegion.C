/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#include "morphingBoxRegion.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

boxBounds morphingBoxRegion::readBounds
(
    const dictionary& dict,
    const scalar defaultValue
)
{
    return
    {
        dict.getOrDefault<scalar>("xMin", -defaultValue),
        dict.getOrDefault<scalar>("xMax", defaultValue),
        dict.getOrDefault<scalar>("yMin", -defaultValue),
        dict.getOrDefault<scalar>("yMax", defaultValue),
        dict.getOrDefault<scalar>("zMin", -defaultValue),
        dict.getOrDefault<scalar>("zMax", defaultValue)
    };
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

morphingBoxRegion::morphingBoxRegion
(
    const dictionary& dict,
    const dictionary* fbPtr
)
{
    // Helper: return dict if it contains key, otherwise fall back to fbPtr
    const auto& src = [&](const word& key) -> const dictionary&
    {
        if (dict.found(key))
        {
            return dict;
        }
        if (fbPtr && fbPtr->found(key))
        {
            return *fbPtr;
        }
        // Return dict and let subsequent lookup issue the error
        return dict;
    };

    // Coordinate system
    cs_ = dict.found("coordinateSystem")
        ? coordinateSystem(dict, "coordinateSystem")
        : (fbPtr && fbPtr->found("coordinateSystem"))
        ? coordinateSystem(*fbPtr, "coordinateSystem")
        : coordinateSystem();

    // Reference frame
    relativeTo_ = src("relativeTo").getOrDefault<word>("relativeTo", "lab");

    // Inner and outer boxes read directly from the region dict
    inner_ = readBounds(src("innerBox").subDict("innerBox"), GREAT);
    outer_ = readBounds(src("outerBox").subDict("outerBox"), VGREAT);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::morphingBoxRegion::updateBlending
(
    const pointField& points,
    scalarField& s
) const
{
    const vector e1(cs_.e1());
    const vector e2(cs_.e2());
    const vector e3(cs_.e3());

    const auto applyRamp = [&]
    (
        const vector& dir,
        const scalar lo,
        const scalar li,
        const scalar ri,
        const scalar ro
    )
    {
        s *= max(min(((points & dir) - lo)/(li - lo), 1.0), 0.0);
        s *= max(min((ro - (points & dir))/(ro - ri), 1.0), 0.0);
    };

    applyRamp(e1, outer_.xMin, inner_.xMin, inner_.xMax, outer_.xMax);
    applyRamp(e2, outer_.yMin, inner_.yMin, inner_.yMax, outer_.yMax);
    applyRamp(e3, outer_.zMin, inner_.zMin, inner_.zMax, outer_.zMax);
}

void Foam::morphingBoxRegion::shiftBoxes(const vector& cor)
{
    if (relativeTo_ != "body") return;

    auto shift = [&](boxBounds& b)
    {
        b.xMin += cor.x(); b.xMax += cor.x();
        b.yMin += cor.y(); b.yMax += cor.y();
        b.zMin += cor.z(); b.zMax += cor.z();
    };

    shift(inner_);
    shift(outer_);
    relativeTo_ = "lab";
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
