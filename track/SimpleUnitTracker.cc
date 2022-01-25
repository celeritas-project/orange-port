//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/SimpleUnitTracker.cc
 * \brief SimpleUnitTracker class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SimpleUnitTracker.hh"

#include <algorithm>
#include <cmath>
#include "base/Assert.hh"
#include "base/Join.hh"
#include "base/Range.hh"
#include "orange/Fuzziness.hh"
#include "orange/surfaces/SurfaceAction.hh"
#include "orange/construct/UnitRegion.hh"
#include "TrackingError.hh"
#include "detail/CellInitializer.hh"
#include "detail/LogicEvaluator.hh"
#include "detail/SurfaceFunctors.hh"
#include "detail/Utils.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// Local helper functions
//---------------------------------------------------------------------------//
namespace
{
//---------------------------------------------------------------------------//
/*!
 * Whether all regions (except possibly the exterior) have the same zorder.
 *
 * The exterior is allowed to have a special value "IMPLICIT_EXTERIOR" if it's
 * non-truncating, indicating that initialization should happen in the closest
 * valid cell. It's not allowed to mask the inside definition, but it's also
 * prevented from initializing during transport.
 *
 * This is essentially the same code as UnitBuilder::is_simple.
 */
bool all_same_zorder(const std::vector<UnitRegion>& regions)
{
    CELER_EXPECT(!regions.empty());
    auto first_interior = regions.begin();
    if (first_interior + 1 != regions.end()
        && first_interior->zorder == ZOrder::implicit_exterior)
    {
        // There are at least two regions: ignore the implicit exterior
        ++first_interior;
    }

    return std::all_of(first_interior + 1,
                       regions.end(),
                       [first_interior](const UnitRegion& reg) {
                           return reg.zorder == first_interior->zorder;
                       });
}

//---------------------------------------------------------------------------//
} // namespace

//---------------------------------------------------------------------------//
/*!
 * Construct from surfaces and regions.
 *
 * All input regions must have equal zorder, since the simple unit tracker
 * doesn't support overlapping regions.
 */
SimpleUnitTracker::SimpleUnitTracker(SurfaceContainer               surfaces,
                                     const std::vector<UnitRegion>& regions)
    : storage_(std::move(surfaces), regions)
{
    CELER_EXPECT(!regions.empty() && all_same_zorder(regions));
}

//---------------------------------------------------------------------------//
/*!
 * Find the local cell and possibly surface ID.
 *
 * TODO: allow initialization on a (single) internal surface? Can we do that
 * only in the case that we're not on a boundary in a higher level (i.e. about
 * to transport out the geometry)?
 */
Initialization SimpleUnitTracker::initialize(LocalState state) const
{
    CELER_EXPECT(state.temp_senses);

    // Find local cells
    auto kd_vols = storage_.find_kdtree_cells(state.pos);

    detail::CellInitializer try_init(storage_.surfaces(), state);

    auto find_in_cell = [&](VolumeId c) -> Initialization {
        Initialization result{};
        if (c == state.cell && state.surface)
        {
            // Cannot cross surface into the same cell
            return result;
        }

        const auto& cell_def = storage_.cells()[c];
        auto        found    = try_init(cell_def);
        if (!found)
        {
            return result;
        }

        result.cell = c;
        if (found.face)
        {
            CELER_ASSERT(found.face < cell_def.faces.size());
            result.surface = cell_def.faces[found.face.unchecked_get()];
            result.sense   = found.sense;
        }

        return result;
    };

    // Loop over possible cells we're inside.
    for (auto kd_volid : kd_vols)
    {
        if (Initialization result = find_in_cell(VolumeId{kd_volid}))
        {
            if (!state.surface && result.surface)
            {
                // Initialized on a boundary but wasn't crossing a surface. Be
                // safe and let the tracking geometry bump it.
                return {};
            }
            return result;
        }
    }

    // Allow missing cells as an unexpected but acceptable error condition:
    // if bumping puts the particle in a valid region then it's not a big
    // deal.
    return {};
}

//---------------------------------------------------------------------------//
/*!
 * Calculate distance-to-intercept for the next surface
 */
Intersection SimpleUnitTracker::intersect(LocalState state) const
{
    CELER_EXPECT(state.cell && state.temp_face_dist && state.temp_senses);
    const auto& cell_def = storage_.cells()[state.cell];

    // Allocate space for all possible surface intersections in the current
    // cell (>= the number of surfaces)
    state.temp_face_dist->resize(cell_def.num_intersections);
    FaceId face = detail::find_face(cell_def.faces, state.surface);
    auto   calc_intersections = make_surface_action(
        storage_.surfaces(), detail::CalcIntersections{state, face});

    // Find all surface intersection distances inside this cell
    for (SurfaceId surface : cell_def.faces)
    {
        CELER_ASSERT(calc_intersections.action().face_dist_iter()
                     != state.temp_face_dist->end());
        calc_intersections(surface);
    }
    CELER_ASSERT(calc_intersections.action().face_dist_iter()
                 == state.temp_face_dist->end());
    CELER_ASSERT(calc_intersections.action().face_idx()
                 == cell_def.faces.size());

    Intersection result;
    if (!(cell_def.flags & detail::CellContainer::INTERNAL_SURFACES))
    {
        // No interior surfaces: closest distance is next boundary
        result = this->simple_intersect(state);
    }
    else
    {
        result = this->complex_intersect(state);
    }

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate normal to the current surface
 */
Real3 SimpleUnitTracker::normal(LocalState state) const
{
    CELER_EXPECT(state.surface);

    return detail::calc_normal(
        storage_.surfaces(), state.pos, state.surface, state.sense);
}

//---------------------------------------------------------------------------//
/*!
 * Access all surfaces
 */
const SurfaceContainer& SimpleUnitTracker::surfaces() const
{
    return storage_.surfaces();
}

//---------------------------------------------------------------------------//
/*!
 * Reconstruct a CSG cell for diagnostic output
 */
UnitRegion SimpleUnitTracker::get_region(VolumeId id) const
{
    return storage_.cells().get_region(id);
}

//---------------------------------------------------------------------------//
/*!
 * Get the face ID of a surface in a cell.
 */
FaceId SimpleUnitTracker::find_face(SurfaceId surf, VolumeId cell) const
{
    CELER_EXPECT(cell && cell.get() < this->num_volumes());
    CELER_EXPECT(surf && surf.get() < this->num_surfaces());
    return detail::find_face(storage_.cells()[cell].faces, surf);
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
/*!
 * Calculate distance to the next boundary for nonreentrant cells
 */
inline Intersection
SimpleUnitTracker::simple_intersect(LocalState state) const
{
    CELER_EXPECT(!state.temp_face_dist->empty());

    // Crossing any surface will leave the cell; perform a linear search for
    // the smallest (but nonzero) distance
    auto crossing = std::min_element(state.temp_face_dist->begin(),
                                     state.temp_face_dist->end(),
                                     detail::CloserFace{});
    CELER_ASSERT(crossing != state.temp_face_dist->end());
    CELER_ASSERT(crossing->second > 0);

    Intersection result;
    if (crossing->second == no_intersection())
        return result;

    // Get surface ID
    const auto& cell_def = storage_.cells()[state.cell];
    result.surface       = cell_def.faces[crossing->first.get()];

    Sense cur_sense{};
    if (result.surface == state.surface)
    {
        // Crossing the same surface that we're currently on (and inside)
        cur_sense = state.sense;
    }
    else
    {
        // XXX find out exactly when this is needed with the new "bump" scheme
        // Calculate current sense with this surface
        SignedSense ss = make_surface_action(
            storage_.surfaces(), detail::CalcSense{state.pos})(result.surface);

        if (NEMESIS_UNLIKELY(ss == SignedSense::on))
        {
            // The same surface may have more than one intersection
            // and we may be positionally "on" this surface while
            // state.surface is unset. This can happen just after moving to
            // the inside of a cylinder/sphere hole: the surface (if any)
            // will be at the lower levels, not at the containing universe.
            Real3 temp_pos = make_vector(state.pos);
            axpy(detail::calc_bump(temp_pos), state.dir, temp_pos);

            ss = make_surface_action(
                storage_.surfaces(),
                detail::CalcSense{temp_pos})(result.surface);
            CELER_ASSERT(ss != SignedSense::on);
        }
        cur_sense = to_sense(ss);
    }

    // Post-surface sense will be on the other side of this cell
    result.sense    = flip_sense(cur_sense);
    result.distance = crossing->second;
    CELER_ENSURE(result && result.distance > 0);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate boundary distance if internal surfaces present.
 *
 * In "complex" cells, crossing a surface can still leave the particle in an
 * "inside" state.
 *
 * We have to iteratively track through all surfaces, in order of minimum
 * distance, to determine whether crossing them will leave us in the cell.
 */
inline Intersection
SimpleUnitTracker::complex_intersect(LocalState state) const
{
    // Get cell definitions for the current cell
    const auto& cell_def = storage_.cells()[state.cell];

    // Find the current face, including possibly coincident faces.
    FaceId face = detail::find_face(cell_def.faces, state.surface);
    detail::calc_senses(
        cell_def, storage_.surfaces(), state.pos, state.temp_senses, &face);
    if (face)
    {
        // On a face: set current sense
        (*state.temp_senses)[face.get()] = state.sense;
    }

    // Current senses should put us inside the cell
    detail::LogicEvaluator is_inside(cell_def.logic);
    CELER_ASSERT(is_inside(*state.temp_senses));

    // Partition into finite nonzero and infinite-or-zero groups.
    auto invalid_iter = std::partition(state.temp_face_dist->begin(),
                                       state.temp_face_dist->end(),
                                       [](const NextFace& fd) {
                                           return fd.second > 0
                                                  && !std::isinf(fd.second);
                                       });

    // Delete the infinites and zeroes
    state.temp_face_dist->erase(invalid_iter, state.temp_face_dist->end());

    // Sort distances (all finite) in ascending order by distance-to-intercept
    std::sort(state.temp_face_dist->begin(),
              state.temp_face_dist->end(),
              [](const NextFace& a, const NextFace& b) {
                  return a.second < b.second;
              });

    // Loop over distances and surface indices to cross.
    // Evaluate the logic expression at each crossing to determine whether
    // we're actually leaving the cell.
    for (const NextFace& crossing : *state.temp_face_dist)
    {
        // Distances should always be positive
        CELER_ASSERT(crossing.second >= 0);

        // Flip sense for the face being crossed
        (*state.temp_senses)[crossing.first.get()].flip();

        if (!is_inside(*state.temp_senses))
        {
            // Flipping this sense puts us outside the cell, so we're done
            // searching.
            Intersection result;
            result.surface = cell_def.faces[crossing.first.get()];
            result.sense = to_sense((*state.temp_senses)[crossing.first.get()]);
            result.distance = crossing.second;
            CELER_ENSURE(result.distance > 0);
            return result;
        }
    }

    // No intersection: perhaps leaving an exterior cell? Possibly geometry
    // error.
    return {};
}

//---------------------------------------------------------------------------//
} // namespace celeritas
