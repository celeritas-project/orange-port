//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/TrackingGeometry.hh
 * \brief TrackingGeometry class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include <vector>
#include "base/FastHashMap.hh"
#include "base/OpaqueIdRange.hh"
#include "orange/Definitions.hh"
#include "TrackingState.hh"
#include "../Definitions.hh"

namespace celeritas
{
class Universe;
//---------------------------------------------------------------------------//
/*!
 * Manage tracking at the highest ("global") level.
 *
 * - Movement within a cell (since no inter-universe traversal is needed)
 * - Position calculation (accounting for within-cell movement)
 * - Direction updates (update global position and propagate to daughters)
 * - Boundary conditions (only meaningful on external surfaces)
 *
 * Within-cell movement is *NOT* propagated into the actual position state:
 * while the particle streams in a straight line, the moved distance is stored
 * independently. This improves efficiency (avoiding vector axpy at every
 * level per movement) and reduces roundoff error (since the position may be
 * orders of magnitude different from the within-cell movement).
 *
 * Sometimes daughter universes will fail to find a position, in which case the
 * tracking geometry will "bump" the internally stored position of the
 * particle. To conserve path length it stores the bump as a *negative*
 * within-cell movement.
 *
 * The expected call sequence is:
 * - \c initialize to set the particle's logical state from a
 *   position/direction
 * - Three functions can be called at any time after initialization:
 *  - \c boundary_state to query whether the particle
 *    is inside or outside or on a reflecting surface
 *  - \c set_direction to change the direction (and invalidate any
 *    distance-to-boundary calculation)
 *  - \c distance_to_boundary to find or return the cached value of the
 *    distance to the next surface
 * - \c intersect (if \c distance_to_boundary hasn't been called) to find the
 *   distance to the next boundary
 * - \c StateView::move_internal to move within the cell
 * - \c cross_surface to move to the boundary
 * - Any of the surface-based operations:
 *   - \c normal to get the normal on the surface
 *   - \c reflect to reflect on an exterior boundary
 * - After moving to a surface, return to \c intersect above.
 *
 * A ray-tracing sequence would call the following sequence of state changes:
 * - \c initialize
 * - \c distance_to_boundary
 * - \c cross_surface
 * - And loop back to intersect if \c boundary_state is still inside.
 */
class TrackingGeometry
{
  public:
    //@{
    //! Public type aliases
    using BoundaryState      = geometria::BoundaryState;
    using SPConstUniverse    = std::shared_ptr<const Universe>;
    using VecUniverse        = std::vector<SPConstUniverse>;
    using MapSurfaceBoundary = FastHashMap<SurfaceId, BoundaryState>;
    using State_t            = TrackingState;
    using UniverseIdRange    = OpaqueIdRange<UniverseId>;
    using size_type          = UniverseId::size_type;
    //@}

    struct Params
    {
        VecUniverse        universes;
        MapSurfaceBoundary boundaries;
    };

  public:
    // Construct with universes and boundary conditions (first is top-level)
    explicit TrackingGeometry(Params params);

    //// ACCESSORS ////

    //! Number of universes
    size_type num_universes() const { return universes_.size(); }

    // Get the universe for this ID
    inline const Universe& get(UniverseId univ_id) const;

    //! Boundary conditions
    const MapSurfaceBoundary& boundaries() const { return boundaries_; }

    //// TRACKER INTERACTION ////

    // Initialize the particle at a position
    void
    initialize(State_t& state, SpanConstReal3 pos, SpanConstReal3 dir) const;

    // Find distance-to-boundary after a surface crossing or direction change
    void intersect(State_t& state) const;

    // Move across the current surface
    void cross_surface(State_t& state) const;

    // Set the direction in absolute coordinates
    void set_direction(State_t& state, SpanConstReal3 dir) const;

    // Calculate the *global* surface normal
    Real3 normal(const State_t& state) const;

    // Reflect on a boundary condition
    void reflect(State_t& state) const;

    // Get or calculate the distance-to-boundary
    inline real_type distance_to_boundary(State_t& state) const;

    //// INTERNAL TO STATE ////

    // Boundary
    BoundaryState boundary_state(const State_t& state) const;

  private:
    //// DATA ////

    // Vector of universes (front is top-level)
    VecUniverse universes_;

    // Boundary conditions
    MapSurfaceBoundary boundaries_;

    //// IMPLEMENTATION ////

    // Initialize or move across boundary
    void initialize_impl(State_t& state) const;

    // Try initializing at a point
    bool try_initialize(State_t& state) const;

    // Whether the state is in the exterior/outside cell
    bool in_exterior(const State_t& state) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "TrackingGeometry.i.hh"
//---------------------------------------------------------------------------//
