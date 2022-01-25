//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/TrackingGeometry.cc
 * \brief TrackingGeometry class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "TrackingGeometry.hh"

#include "base/FixedViewArray.hh"
#include "base/Range.hh"
#include "base/VectorFunctions.hh"
#include "orange/Fuzziness.hh"
#include "StateView.hh"
#include "TrackingError.hh"
#include "Universe.hh"
#include "detail/LocalStateBuilder.hh"
#include "detail/Transform.hh"

using make_fixed_view;

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with universes (first is top-level)
 */
TrackingGeometry::TrackingGeometry(Params params)
    : universes_(std::move(params.universes))
    , boundaries_(std::move(params.boundaries))
{
    CELER_EXPECT(!universes_.empty());
#ifdef REQUIRE_ON
    // Check that universe IDs match index position
    for (auto i : range(universes_.size()))
    {
        CELER_EXPECT(universes_[i]);
        CELER_EXPECT(universes_[i]->id().get() == i);
    }
    // Check boundary conditions
    SurfaceId::size_type top_num_surfaces = universes_.front()->num_surfaces();
    for (const auto& surfid_boundary : boundaries_)
    {
        CELER_EXPECT(surfid_boundary.first.get() < top_num_surfaces);
        CELER_EXPECT(surfid_boundary.second == geometria::REFLECT);
    }
#endif
    // Outermost universe must have both interior and exterior
    CELER_EXPECT(universes_.front()->num_volumes() >= 2);
}

//---------------------------------------------------------------------------//
// TRACKER INTERACTION
//---------------------------------------------------------------------------//
/*!
 * Initialize the particle at a position.
 *
 * To aid in user error conditions, the state must always kept up-to-date with
 * the current level/coordinates being tested.
 */
void TrackingGeometry::initialize(State_t&       state,
                                  SpanConstReal3 pos,
                                  SpanConstReal3 dir) const
{
    // Initialize data using state view
    StateView state_view(state);
    state_view = {pos, dir};

    // Initialize tracking
    this->initialize_impl(state);
}

//---------------------------------------------------------------------------//
/*!
 * Find the next surface along the current direction.
 *
 * TODO: discriminate between "unchecked" intersections and "known infinite"
 * intersections so we don't have to repeatedly check infinite distances.
 */
void TrackingGeometry::intersect(State_t& state) const
{
    CELER_EXPECT(state);
    CELER_EXPECT(state.movement <= 0);
    CELER_EXPECT(!StateView(state).has_next_distance());

    const real_type bump_rel = 1.0 + fuzziness().bump_rel();

    // Initialize 'next level' to top universe, so that an infinite distance
    // at every level will be a movement on the top level universe
    int       next_level    = 0;
    real_type next_distance = no_intersection();

    detail::LocalStateBuilder build_local_state(&state);

    // Local state should *not* have surface/sense except at the lowest level
    for (int level : range(state.level + 1))
    {
        Intersection& next = state.level_state[level].next;
        if (!next)
        {
            // Level needs to be recalculated
            UniverseId univ_id = state.level_state[level].univ_id;
            CELER_ASSERT(univ_id < universes_.size());
            next = universes_[univ_id.get()]->intersect(
                build_local_state(level));
            CELER_ASSERT(!next || next.distance > 0);
        }

        // Use *soft less than* to prefer staying in higher level universe on
        // coincident boundaries: effectively make deeper-level distance look
        // "a little more infinite" so that its outer boundary is less
        // preferable to the parent level.
        if (next.distance * bump_rel < next_distance)
        {
            next_level    = level;
            next_distance = next.distance;
        }
    }

    state.next_level = next_level;
    CELER_ENSURE(StateView(state).has_next_distance());
}

//---------------------------------------------------------------------------//
/*!
 * Move across the current surface
 */
void TrackingGeometry::cross_surface(State_t& state) const
{
    CELER_EXPECT(StateView(state).has_next_distance());
    // New level is "next" level: either at the highest level (move to adjacent
    // cell in local universe, possibly pushing into another level) *or* at a
    // higher level (popping out of the universe).
    state.level      = state.next_level;
    state.next_level = -1;

    // Update position so the local state is at the intersection point
    const auto next = state.level_state[state.level].next;
    if (next.distance > 0)
    {
        state.movement = next.distance;
        StateView(state).apply_movement();
    }

    state.on_surface.surface = next.surface;
    state.on_surface.sense   = next.sense;
    state.on_surface.level   = state.level;

    // Move to next cell, dive into daughters as necessary, and clear "next".
    this->initialize_impl(state);
    CELER_ENSURE(!StateView(state).has_next_distance());
}

//---------------------------------------------------------------------------//
/*!
 * Set the direction in absolute coordinates
 */
void TrackingGeometry::set_direction(State_t& state, SpanConstReal3 dir) const
{
    if (state.movement > 0)
    {
        // Update point along current direction
        StateView(state).apply_movement();
    }
    else if (state.movement < 0)
    {
        // We're change direction before moving along a bump to a known-good
        // position. Roll the bump distance into the overall correction term
        // to ensure that our logical position is consistent (so we don't
        // accidentally move backward out of our cell).
        axpy(state.movement, state.level_state[0].dir, state.correction);
        state.movement   = 0;
        state.on_surface = {};
    }

    // Loop over all levels, setting direction in current level, rotating into
    // daughter's reference frame, and getting pointer to the daughter
    // universe.
    Real3           local_dir  = make_vector(dir);
    const Universe* local_univ = universes_.front().get();
    for (int level : range(state.level))
    {
        state.level_state[level].dir = local_dir;

        const auto* daughter
            = local_univ->daughter(state.level_state[level].cell);
        CELER_ASSERT(daughter);
        local_univ = daughter->universe.get();

        // Apply rotation to update local_dir to daughter level
        if (daughter->transform.has_rotation())
        {
            daughter->transform.rotate_to_daughter(local_dir);
        }

        // Clear distance-to-boundary on all levels
        state.level_state[level].next = {};
    }

    // Update deepest level direction and intersection
    CELER_ASSERT(universes_[state.level_state[state.level].univ_id.get()].get()
                 == local_univ);
    state.level_state[state.level].dir  = local_dir;
    state.level_state[state.level].next = {};

    state.next_level = -1;

    CELER_ENSURE(state.movement == 0);
    CELER_ENSURE(!StateView(state).has_next_distance());
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the surface normal in the global frame
 */
Real3 TrackingGeometry::normal(const State_t& state) const
{
    CELER_EXPECT(state.on_surface);
    CELER_EXPECT(state.on_surface.level >= 0
                 && state.on_surface.level <= state.level);

    int   level = state.on_surface.level;
    Real3 result_dir;

    // Calculate local normal from universe with surface
    {
        auto       local_state = detail::LocalStateBuilder(&state)(level);
        UniverseId univ        = state.level_state[level].univ_id;

        result_dir = this->get(univ).normal(local_state);
        --level;
    }

    // Propagate normal direction upward through levels: first iteration is
    // above, last iteration is for level zero.
    for (; level >= 0; --level)
    {
        const auto& level_state = state.level_state[level];
        UniverseId  univ        = level_state.univ_id;

        const auto* daughter = this->get(univ).daughter(level_state.cell);
        CELER_ASSERT(daughter);
        if (daughter->transform.has_rotation())
        {
            daughter->transform.rotate_to_parent(result_dir);
        }
    }
    return result_dir;
}

//---------------------------------------------------------------------------//
// INTERNAL TO STATE
//---------------------------------------------------------------------------//
/*!
 * Determine the boundary state.
 *
 * Any point inside the exterior boundary is inside. When the particle is in
 * the outside cell but on the surface, it might have a special boundary. If
 * outside but not on a surface, it is always outside.
 */
auto TrackingGeometry::boundary_state(const State_t& state) const
    -> BoundaryState
{
    CELER_EXPECT(state);
    if (state.level > 0)
    {
        // Not the top-level universe; therefore inside the geometry
        return geometria::inside;
    }
    if (this->in_exterior(state))
    {
        if (state.on_surface && state.movement <= 0)
        {
            // On a top-level surface: see if it's a special boundary
            auto iter = boundaries_.find(state.on_surface.surface);
            if (iter != boundaries_.end())
            {
                // Return special boundary condition
                return iter->second;
            }
        }
        return geometria::outside;
    }
    return geometria::inside;
}

//---------------------------------------------------------------------------//
/*!
 * Reflect on a boundary.
 *
 * This may only be called if the particle is on the exterior of a
 * reflecting boundary. It's like a cross-surface but with a distance of zero
 * on the current face rather than a "next" face.
 *
 * Reflection: \f[
  \Omega' = \Omega - 2 (\Omega \cdot n) n
  \f]
 * The sign of \em n doesn't matter because it cancels out, so it's not
 * important whether the low-level tracker provides an outward- or
 * inward-facing normal.
 */
void TrackingGeometry::reflect(State_t& state) const
{
    CELER_EXPECT(state.level == 0);
    CELER_EXPECT(state.movement <= 0);
    CELER_EXPECT(state.on_surface);

    if (NEMESIS_UNLIKELY(state.movement < 0))
    {
        // Had to bump to cross this boundary. Roll the bump distance into the
        // overall correction term to ensure that our logical position is
        // consistent (so we don't accidentally move backward out of our cell).
        axpy(state.movement, state.level_state[0].dir, state.correction);
        state.movement = 0;
    }

    // Reflect the direction across the surface normal:
    // dir' = -2 (dir . normal) * normal + dir
    auto normal = this->normal(state);
    auto dir    = make_fixed_view(state.level_state[0].dir);
    axpy(-2 * dot_product(normal, dir), normal, dir);

    // Next sense is flipped because we're staying on the same side of the
    // surface instead of crossing it.
    state.on_surface.sense = flip_sense(state.on_surface.sense);

    // Distance-to-boundary shouldn't have been called..
    CELER_ASSERT(!state.level_state[0].next);

    // Leave the outside cell, cross back inside (using the current
    // cell/surface/sense information)
    return this->initialize_impl(state);
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * Set local tracking information and move to daughters.
 *
 * This is a very slight generalization of initialize + cross_surface.
 */
void TrackingGeometry::initialize_impl(State_t& state) const
{
    CELER_EXPECT(state);

    /*!
     * Typical initialization.
     *
     * Try actual position and stated level, descending downward through
     * universes until the final one is found.
     */
    if (this->try_initialize(state))
    {
        // Usual case: initialization succeeded
        return;
    }

    Real3     orig_pos = state.level_state[0].pos;
    real_type dx       = detail::calc_bump(orig_pos);

    /*!
     * Along-direction bump.
     *
     * Move along the particle's direction by the bump distance and see if
     * initialization succeeds. This is preferable to the "corner" bumping
     * below because the internal movement logic is able to seamlessly
     * correct for the bump distance at every step or direction change.
     *
     * TODO: add "is inside cell" logic to tracker so that we can loop upward
     * from lowest level rather than reinitializing from scratch every time
     *
     * TODO: if bump succeeded, loop backward through levels, looking for
     * nearest surface in opposite direction so that we keep an accurate
     * "universe"/"on surface"
     */
    state.level      = 0;
    state.next_level = -1;
    state.on_surface = {};

    axpy(dx,
         make_fixed_view(state.level_state[0].dir),
         make_fixed_view(state.level_state[0].pos));
    if (this->try_initialize(state))
    {
        state.movement -= dx;
        return;
    }

    /*!
     * Nearby-position bump.
     *
     * This is most commonly needed when initializing exactly on *and along* a
     * planar surface. Start in the positive directions first just because
     * users and unit test writers are more likely to raytrace along the
     * positive axes (so at lower corners this tends to put us "inside" the
     * geometry instead of "outside").
     */
    // Move at most a radial distance of the bump distance
    constexpr real_type inv_sqrt_three = 0.5773502691896257;
    dx *= inv_sqrt_three;

    Real3 delta;
    for (int corner_bits : range(1 << 3))
    {
        for (int ax : range(3))
        {
            delta[ax] = (corner_bits & (1 << ax) ? -dx : dx);
        }

        // Calculate bump along the direction of flight
        real_type along_bump = dot_product(delta, state.level_state[0].dir);
        if (along_bump < 0)
        {
            // Never bump backward!!
            continue;
        }

        // Try initializing from scratch at this position
        state.level      = 0;
        state.next_level = -1;
        state.on_surface = {};
        for (int ax : range(3))
        {
            state.level_state[0].pos[ax] = orig_pos[ax] + delta[ax];
        }

        if (this->try_initialize(state))
        {
            // Initialization succeeded: split movement into "along-direction"
            // movement and azimuthal correction
            state.movement -= along_bump;

            // Subtract along-particle component from correction
            for (int ax : range(3))
            {
                state.correction[ax]
                    += -delta[ax] + state.level_state[0].dir[ax] * along_bump;
            }

            return;
        }
    }

    // TODO: better job of restoring state for error reporting
    state.level_state[0].pos = orig_pos;
    ORANGE_TRACKING_ASSERT(false, NoCellError());
}

bool TrackingGeometry::try_initialize(State_t& state) const
{
    CELER_EXPECT(state.level >= 0);
    detail::LocalStateBuilder build_local_state(&state);

    const Universe* universe = nullptr;

    // Initialize starting level
    {
        auto& level_state = state.level_state[state.level];
        universe          = &this->get(level_state.univ_id);
        level_state.next  = {};
    }

    while (true)
    {
        CELER_ASSERT(universe);
        LocalState parent = build_local_state(state.level);
        Initialization     init   = universe->initialize(parent);
        if (!init)
        {
            // Initialization failed! Revert to the higher level (mainly
            // important for debugging), then try bumping and recursing again.
            --state.level;
            return false;
        }

        // Save successful initialization for this level
        CELER_ASSERT(init.cell);
        state.level_state[state.level].cell = init.cell;

        if (state.level == state.on_surface.level)
        {
            // Update surface initialization at this level
            state.on_surface.surface = init.surface;
            state.on_surface.sense   = init.sense;
        }
        else
        {
            // With the current implementation we should never end up on a
            // surface unless we're known to be crossing it from a previous
            // "intersect".
            CELER_ASSERT(!init.surface);
        }

        if (const auto* daughter = universe->daughter(init.cell))
        {
            // Traverse into the daughter universe
            ++state.level;
            CELER_ASSERT(state.level < state.max_level);
            auto& level_state = state.level_state[state.level];
            universe          = daughter->universe.get();
            // Set pos/dir in new level from transformed parent coords
            detail::transform_to_daughter(parent.pos,
                                          parent.dir,
                                          daughter->transform,
                                          make_fixed_view(level_state.pos),
                                          make_fixed_view(level_state.dir));

            // Set/clear properties for new level
            level_state.cell    = {};
            level_state.univ_id = universe->id();
            level_state.next    = {};
        }
        else
        {
            // "Leaf" in the universe tree: found a physical material, so no
            // more descending into daughters.
            break;
        }
    }

    return true;
    CELER_ENSURE(!StateView(state).has_next_distance());
}

//---------------------------------------------------------------------------//
/*!
 * Whether the state is in the exterior/outside cell.
 *
 * This is true if and only if the particle is in the outermost universe and
 * inside the lowest-numbered cells. The "outside" cell is zero rather than
 * (num_volumes - 1), because masked regions are ordered in *decreasing*
 * priority and the "EXTERIOR" zorder is higher than anything.
 */
bool TrackingGeometry::in_exterior(const State_t& state) const
{
    return (state.level == 0) && (state.level_state[0].cell == VolumeId{0});
}

//---------------------------------------------------------------------------//
} // namespace celeritas
