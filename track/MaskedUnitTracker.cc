//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/MaskedUnitTracker.cc
 * \brief MaskedUnitTracker class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "MaskedUnitTracker.hh"

#include <algorithm>
#include <functional>
#include "base/FixedViewArray.hh"
#include "orange/Fuzziness.hh"
#include "orange/surfaces/SurfaceAction.hh"
#include "orange/construct/UnitRegion.hh"
#include "TrackingError.hh"
#include "detail/CellInitializer.hh"
#include "detail/FoundCells.hh"
#include "detail/LogicEvaluator.hh"
#include "detail/SurfaceFunctors.hh"
#include "detail/Utils.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct from surfaces and regions.
 */
MaskedUnitTracker::MaskedUnitTracker(SurfaceContainer               surfaces,
                                     const std::vector<UnitRegion>& regions)
    : storage_(std::move(surfaces), regions), surface_groups_(storage_)
{
}

//---------------------------------------------------------------------------//
/*!
 * Find the local cell and possibly surface ID.
 */
Initialization MaskedUnitTracker::initialize(LocalState state) const
{
    // Find local cells
    auto kd_vols = storage_.find_kdtree_cells(state.pos);

    detail::FoundCells      found_cells;
    detail::CellInitializer try_init(storage_.surfaces(), state);

    // Loop over possible cells we're inside. Cells are ordered by decreasing Z
    // number.
    for (auto kd_volid : kd_vols)
    {
        VolumeId test_cell{kd_volid};
        if (state.surface && test_cell == state.cell)
        {
            // Cannot cross surface into the same cell
            continue;
        }

        const auto& cell_def = storage_.cells()[test_cell];

        if (!found_cells.empty() && cell_def.zorder != found_cells.zorder())
        {
            // We already found a cell on a surface at a higher zorder, so
            // there are no other potential overlapping cells.
            CELER_ASSERT(cell_def.zorder < found_cells.zorder());
            break;
        }

        if (auto found = try_init(cell_def))
        {
            // Cell is found! Save the cell ID and zorder
            found_cells.insert(test_cell, cell_def.zorder, found);
        }
    }

    if (state.cell && found_cells.size() == 1
        && state.cell != found_cells.cell_view().front())
    {
        // Exited an old cell into a new cell. If the Z order is the same,
        // make sure the new position is *not* inside the old cell.
        const auto& cell_def = storage_.cells()[state.cell];
        if (cell_def.zorder == found_cells.zorder() && try_init(cell_def))
        {
            found_cells.insert(state.cell, cell_def.zorder, FoundFace{true});
        }
    }

    ORANGE_TRACKING_ASSERT(found_cells.size() <= 1,
                           OverlappingCellError{found_cells.cell_view()});

    if (found_cells.empty())
    {
        // Allow missing cells as an unexpected but acceptable error condition:
        // if bumping puts the particle in a valid region then it's not a big
        // deal.
        return {};
    }

    Initialization result;
    result.cell = found_cells.cell_view().front();
    auto face   = found_cells.face();
    if (face)
    {
        if (!state.surface)
        {
            // Initialized on a boundary but wasn't crossing a surface. Be safe
            // and let the tracking geometry bump it.
            return {};
        }
        const auto& cell_def = storage_.cells()[result.cell];
        result.surface       = cell_def.faces[face.unchecked_get()];
        result.sense         = found_cells.sense();
    }

    // Cell should change if crossing a surface
    CELER_ENSURE((result.cell != state.cell) || !state.surface);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate distance-to-intercept for the next surface.
 *
 * Currently we only find the distance to the nearest cell with order greater
 * than or equal to the current cell. Collision detection should be handled
 * automatically in the `initialize` when crossing the surface.
 *
 * Possible intersections are:
 * - Any surface in a higher zorder cell (we must check whether the cell
 *   on the other side is higher-zorder when crossing)
 * - A surface connected to the current cell (leaving it could put us in a
 *   lower zorder cell)
 *
 * TODO: calculate all intersections and all senses once, then reconstruct
 * per-cell and per-face from those lists. Otherwise we're recalculating the
 * same surface multiple times (if it's present in multiple z orders).
 */
Intersection MaskedUnitTracker::intersect(LocalState state) const
{
    CELER_EXPECT(state.cell);

    const zorder_int cur_zorder = storage_.cells()[state.cell].zorder;
    const real_type  bump_dist  = detail::calc_bump(state.pos);
    Intersection     result;

    // Loop over zorders (in decreasing order) in this unit
    for (const auto& surf_group : surface_groups_)
    {
        if (surf_group.zorder < cur_zorder)
        {
            // We're past potential surfaces that could apply
            break;
        }
        // Find current surface in the list of similarly-zorder surfaces
        FaceId cur_face = detail::find_face(surf_group.faces, state.surface);

        // Allocate space for all possible surface intersections
        state.temp_face_dist->resize(surf_group.num_intersections);
        auto calc_intersections = make_surface_action(
            storage_.surfaces(), detail::CalcIntersections{state, cur_face});

        // Find all surface intersection distances inside this group
        for (SurfaceId sid : surf_group.faces)
        {
            CELER_ASSERT(calc_intersections.action().face_dist_iter()
                         != state.temp_face_dist->end());
            calc_intersections(sid);
        }
        CELER_ASSERT(calc_intersections.action().face_dist_iter()
                     == state.temp_face_dist->end());
        CELER_ASSERT(calc_intersections.action().face_idx()
                     == surf_group.faces.size());

        // Loop over surfaces closer than the current
        for (const NextFace& crossing : *state.temp_face_dist)
        {
            if (crossing.second <= 0)
            {
                // Skip coincident or negative surfaces
                continue;
            }
            else if (crossing.second >= result.distance)
            {
                // Skip this surface crossing because it's further away than an
                // already found possibility
                continue;
            }

            // Calculate position just past the surface
            Real3 pos = make_vector(state.pos);
            axpy(crossing.second + bump_dist, state.dir, make_fixed_view(pos));

            // Loop over cells connected to this surface that have the surface
            // group's zorder
            SurfaceId sid = surf_group.faces[crossing.first.get()];
            for (VolumeId::size_type c : storage_.get_connected_cells(sid))
            {
                const auto& cell_def = storage_.cells()[VolumeId{c}];
                if (cell_def.zorder != surf_group.zorder)
                {
                    // Only testing cells of the current zorder surface group
                    continue;
                }

                // Test senses on the new cell: either inside the tested cell
                // or leaving the current cell. Assume that bumped past this
                // target surface we are not on a face (or if we are, it
                // doesn't matter) so set face to nullptr.
                detail::calc_senses(cell_def,
                                    storage_.surfaces(),
                                    pos,
                                    state.temp_senses,
                                    nullptr);

                detail::LogicEvaluator is_inside(cell_def.logic);
                if (is_inside(*state.temp_senses)
                    != (VolumeId{c} == state.cell))
                {
                    // Either we are *no longer* in the starting cell or *are*
                    // in a new cell.
                    result.distance = crossing.second;
                    result.surface  = sid;

                    // Calculate the sense of the crossed face in the new cell
                    auto on_face = detail::find_face(cell_def.faces, sid);
                    CELER_ASSERT(on_face < state.temp_senses->size());
                    result.sense
                        = to_sense((*state.temp_senses)[on_face.get()]);

                    break;
                }
            }
        }
    }

    CELER_ENSURE(!result || result.distance > 0);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate normal to the current surface
 */
Real3 MaskedUnitTracker::normal(LocalState state) const
{
    CELER_EXPECT(state.surface);

    return detail::calc_normal(
        storage_.surfaces(), state.pos, state.surface, state.sense);
}

//---------------------------------------------------------------------------//
/*!
 * Access all surfaces
 */
const SurfaceContainer& MaskedUnitTracker::surfaces() const
{
    return storage_.surfaces();
}

//---------------------------------------------------------------------------//
/*!
 * Reconstruct a CSG cell for diagnostic output
 */
UnitRegion MaskedUnitTracker::get_region(VolumeId id) const
{
    return storage_.cells().get_region(id);
}

//---------------------------------------------------------------------------//
/*!
 * Get the face ID of a surface in a cell.
 */
FaceId MaskedUnitTracker::find_face(SurfaceId surf, VolumeId cell) const
{
    CELER_EXPECT(cell && cell.get() < this->num_volumes());
    CELER_EXPECT(surf && surf.get() < this->num_surfaces());
    auto face = detail::find_face(storage_.cells()[cell].faces, surf);

    return face;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
