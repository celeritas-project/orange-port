//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/RectArrayTracker.cc
 * \brief RectArrayTracker class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RectArrayTracker.hh"

#include <type_traits>
#include "base/Face.hh"
#include "base/Range.hh"
#include "detail/GridFinder.hh"

using RectFace = celeritas::RectArrayTracker::Face;

namespace celeritas
{
namespace
{
//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
//! Convert a surface ID to a cuboid face
inline RectFace to_face(SurfaceId s)
{
    CELER_EXPECT(s < RectFace::num_faces());
    return RectFace(s.unchecked_get());
}

//---------------------------------------------------------------------------//
//! Convert a cuboid face to a surface ID
inline SurfaceId to_surface(RectFace s)
{
    CELER_EXPECT(s);
    return SurfaceId(s.unchecked_get());
}

//---------------------------------------------------------------------------//
} // namespace

//---------------------------------------------------------------------------//
/*!
 * Construct a rect array tracker from a cartesiain grid.
 */
RectArrayTracker::RectArrayTracker(Grid_t grid) : grid_(std::move(grid))
{
    CELER_EXPECT(grid_);
}

//---------------------------------------------------------------------------//
/*!
 * Find the new local cell.
 *
 * TODO: with added "sense" initialize/intersect we could probably simplify
 * this. Crossing initialization will have 'outside' sense, other
 * initialization will be 'inside'.
 *
 * TODO: if we initialize outside of the "valid" region, return
 * Initialization{} so that the higher level geometry can bump and try
 * somewhere else.
 */
Initialization RectArrayTracker::initialize(LocalState state) const
{
    if (state.surface)
    {
        CELER_EXPECT(state.cell < grid_.num_volumes());
        return this->cross_surface(this->cell_to_coords(state.cell),
                                   to_face(state.surface));
    }
    else
    {
        return this->initialize_interior(state);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Calculate distance-to-intercept for the next surface.
 */
Intersection RectArrayTracker::intersect(LocalState state) const
{
    CELER_EXPECT(state.cell);

    Grid_t::DimVector coords = this->cell_to_coords(state.cell);
    Face cur_face            = state.surface ? to_face(state.surface) : Face{};

    Intersection result;
    result.sense = Sense::outside;
    for (auto ax : range<FaceId::size_type>(3))
    {
        // Ignore any stationary axis
        if (state.dir[ax] == 0.0)
            continue;

        // Positive direction -> leaving negative surface, will intersect
        // positive
        Face target_face(ax, state.dir[ax] > 0);
        if (cur_face == target_face)
        {
            // Don't intersect surface being crossed
            continue;
        }

        // Find the target edge space position
        int target_coord = int(coords[ax]) + target_face.is_positive();

        // Calculate intersection distance
        const auto& edges = grid_.edges(ax);
        CELER_ASSERT(target_coord >= 0 && target_coord < edges.size());
        real_type dist = (edges[target_coord] - state.pos[ax]) / state.dir[ax];
        if (dist > 0 && dist < result.distance)
        {
            result.distance = dist;
            result.surface  = to_surface(target_face);
        }
    }

    CELER_ENSURE(!result || result.distance > 0);
    CELER_ENSURE(!result.surface || result.surface.get() < 6);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate normal to the current surface
 */
Real3 RectArrayTracker::normal(LocalState state) const
{
    CELER_EXPECT(state.surface);

    // Convert face ID into a cartesian face indexer
    auto face = to_face(state.surface);
    CELER_ASSERT(face && face.axis() < 3
                 && (face.normal() == 1 || face.normal() == -1));

    Real3 result{0, 0, 0};
    result[face.axis()] = face.normal();
    return result;
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Convert cell ID to ijk coords
 */
auto RectArrayTracker::cell_to_coords(VolumeId cell) const -> DimVector
{
    CELER_EXPECT(cell < grid_.num_volumes());
    return grid_.index(cell.unchecked_get());
}

//---------------------------------------------------------------------------//
/*!
 * Convert IJK to cell ID
 */
auto RectArrayTracker::coord_to_cell(const DimVector& ijk) const -> VolumeId
{
    CELER_EXPECT(ijk.all_lt(grid_.num_cells_dims()));
    return VolumeId(grid_.index(ijk));
}

//---------------------------------------------------------------------------//
/*!
 * Initialize inside the array, not previously on a surface
 */
Initialization
RectArrayTracker::initialize_interior(LocalState state) const
{
    CELER_EXPECT(!state.surface);

    // Find coords from spatial location.
    DimVector coords;
    Face      on_face;
    for (auto ax : range(def::END_IJK))
    {
        detail::GridFinder<real_type> find_cell(make_span(grid_.edges(ax)));

        auto found = find_cell(state.pos[ax]);
        coords[ax] = found.cell;
        if (found.edge)
        {
            // Reject initialization on edges
            return {};
        }
    }

    Initialization result;
    result.cell = this->coord_to_cell(coords);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Initialize in the cell across the current surface.
 */
Initialization
RectArrayTracker::cross_surface(DimVector coords, Face face) const
{
    CELER_EXPECT(face);

    // Use face axis/positivity to move to adjacent cell, clamping to
    // nearest interior cell.  Use *signed* arithmetic to find new coordinate
    // to avoid overflow
    int next_coord = int(coords[face.axis()]) + face.normal();
    if (0 <= next_coord && next_coord < int(grid_.num_cells_along(face.axis())))
    {
        // Still inside the grid: update coordinate and flip to the
        // opposing face
        coords[face.axis()] = next_coord;
        face                = face.opposite();
    }

    Initialization result;
    result.cell    = this->coord_to_cell(coords);
    result.surface = to_surface(face);
    result.sense   = Sense::inside;

    return result;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
