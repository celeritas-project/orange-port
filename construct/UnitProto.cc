//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/UnitProto.cc
 * \brief UnitProto class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "UnitProto.hh"

#include <algorithm>
#include "base/Future.hh"
#include "orange/track/SimpleUnitTracker.hh"
#include "orange/track/MaskedUnitTracker.hh"
#include "IntersectionShape.hh"
#include "UnitBuilder.hh"

namespace
{
struct is_true
{
    template<class T>
    bool operator()(const T& obj)
    {
        return static_cast<bool>(obj);
    }
};
} // namespace

namespace celeritas
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * Construct a shape representing the inside of a "well-connected" hole.
 *
 * This is to be used for KENO-V.a geometries as well as Geometria-style holes.
 * The resulting shapes should be excluded from appropriate regions.
 */
auto UnitProto::make_hole_shape(const Hole& hole) -> SPConstShape
{
    CELER_EXPECT(hole.proto);

    return shape_from_rdv(hole.proto->interior(),
                          hole.transform,
                          hole.md ? hole.md : hole.proto->metadata());
}

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Constructor
 */
UnitProto::UnitProto(Params params) : data_(std::move(params))
{
    Insist(data_.md, "Missing unit metadata");
    Insist(data_.boundary,
           "Incomplete interior definition for unit '" << data_.md.name()
                                                       << "'");
    Insist(!data_.media.empty() || !data_.arrays.empty()
               || !data_.holes.empty(),
           "No media/arrays/holes are defined for unit '" << data_.md.name());
    Insist(std::all_of(data_.media.begin(), data_.media.end(), is_true{}),
           "Incomplete media definitions for unit '" << data_.md.name() << "'");
    Insist(std::all_of(data_.holes.begin(), data_.holes.end(), is_true{}),
           "Incomplete hole definitions for unit '" << data_.md.name() << "'");
    Insist(std::is_sorted(data_.holes.begin(),
                          data_.holes.end(),
                          [](const Hole& lhs, const Hole& rhs) {
                              return lhs.zorder > rhs.zorder;
                          }),
           "Holes must be provided in decreasing z order");
    CELER_ASSERT(data_.holes.empty()
                 || data_.holes.back().zorder <= data_.holes.front().zorder);
}

//---------------------------------------------------------------------------//
//! Default destructor
UnitProto::~UnitProto() = default;

//---------------------------------------------------------------------------//
/*!
 * Construct the tracker, metadata, etc.
 *
 * Regarding boundary conditions. "Params" is for how the unit is *defined*:
 *
 *  ========== ==========================================================
 *   Params     Description
 *  ========== ==========================================================
 *   Implicit   Boundary implicitly truncates interior (KENO)
 *   Explicit   Interior is explicitly declared inside boundary (Oberon)
 *  ========== ==========================================================
 *
 * "Build args" is for the context in which the unit is *built*:
 *
 *  ========== ==========================================================
 *   Args       Description
 *  ========== ==========================================================
 *   Implicit   Boundary is already truncated by higher-level unit
 *   Explicit   Boundary must explicitly be represented (global unit)
 *  ========== ==========================================================
 *
 * Args.implicit_boundary is false if and only if the unit is "global".
 *
 *  === === ==========================================================
 *   A   P   Result
 *  === === ==========================================================
 *   I   I  IMPLICIT_EXTERIOR: Higher-level universe truncates
 *   I   X  IMPLICIT_EXTERIOR: Higher-level universe truncates
 *   X   I  EXTERIOR: Global unit that truncates other regions
 *   X   X  MEDIA: Global unit with well-connected exterior
 *  === === ==========================================================
 *
 *  When "args" is implicit, the boundary shape's CSG subtree expression can be
 *  replaced with TRUE wherever it occurs inside the shapes.
 */
auto UnitProto::build(BuildArgs args) const -> BuildResult
{
    BuildResult result;

    // Construct
    UnitBuilder build_unit;
    build_unit.reserve(data_.media.size() + data_.arrays.size()
                       + data_.holes.size());

    // Build exterior cell
    build_unit.exterior(data_.boundary.interior,
                        args.implicit_boundary ? ZOrder::implicit_exterior
                        : data_.boundary.implicit_boundary ? ZOrder::exterior
                                                           : ZOrder::media,
                        data_.boundary.md);

    // Build holes
    for (const Hole& hole : data_.holes)
    {
        CELER_ASSERT(hole.zorder >= ZOrder::media);
        // Construct shape
        SPConstShape hole_shape = UnitProto::make_hole_shape(hole);
        // Add region
        auto volume_id
            = build_unit.region({{neg, hole_shape}}, hole.zorder, hole.md);
        // Add daughter universe
        result.daughters.push_back(
            {volume_id, Daughter{hole.proto, hole.transform}});
    }

    // Build arrays
    for (const Array& arr : data_.arrays)
    {
        // Add region
        auto volume_id = build_unit.region(arr.interior, ZOrder::media, arr.md);
        // Add daughter universe
        result.daughters.push_back(
            {volume_id, Daughter{arr.proto, arr.transform}});
    }

    // Build media
    for (const Media& med : data_.media)
    {
        // Add region
        auto volume_id = build_unit.region(
            med.interior, ZOrder::media, med.md, med.volume);
        // Add matid
        result.matids.push_back({volume_id, med.matid});
    }

    // Construct the tracker and metadata
    auto unit_components = build_unit(data_.md);
    if (unit_components.md->is_simple() && !data_.otf_error_checking)
    {
        result.tracker = make_unique<SimpleUnitTracker>(
            std::move(unit_components.surfaces), unit_components.regions);
    }
    else
    {
        result.tracker = make_unique<MaskedUnitTracker>(
            std::move(unit_components.surfaces), unit_components.regions);
    }
    result.md = std::move(unit_components.md);

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Get the boundary definition for defining a hole in a higher level
 */
auto UnitProto::interior() const -> const RegionVec&
{
    return data_.boundary.interior;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
