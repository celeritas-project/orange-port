//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/StateView.i.hh
 * \brief StateView inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Assert.hh"
#include "base/Range.hh"
#include "base/VectorFunctions.hh"
#include "detail/Utils.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTION
//---------------------------------------------------------------------------//
/*!
 * Constructor takes reference to the state.
 *
 * To facilitate integration with the rest of the CPU code, which passes around
 * mutable/const pointers to a single object, the constructor accepts the state
 * data and const-casts it. The const/mutable semantics for this class
 * member functions determine whether the data is actually modified, so if the
 * caller constructs a \c StateView with a const state, it should not call any
 * of the mutating methods.
 */
StateView::StateView(const TrackingState& state)
    : state_(const_cast<TrackingState&>(state))
{
}

//---------------------------------------------------------------------------//
// INITIALIZATION
//---------------------------------------------------------------------------//
/*!
 * Initialize the state
 */
StateView& StateView::operator=(const Initializer_t& init)
{
    CELER_EXPECT(soft_unit_vector(init.dir));

    state_.level                  = 0;
    state_.next_level             = -1;
    state_.level_state[0].pos     = make_vector(init.pos);
    state_.level_state[0].dir     = make_vector(init.dir);
    state_.level_state[0].cell    = {};
    state_.level_state[0].univ_id = UniverseId{0};
    state_.level_state[0].next    = {};
    state_.on_surface             = {};
    state_.movement               = 0;
    state_.correction             = {0, 0, 0};

    CELER_ENSURE(!this->has_next_distance());
    return *this;
}

//---------------------------------------------------------------------------//
// ACCESSORS
//---------------------------------------------------------------------------//
/*!
 * Get all the state for the highest (global) level.
 */
auto StateView::global_state() const -> const LevelState&
{
    CELER_EXPECT(state_);
    return state_.level_state[0];
}

//---------------------------------------------------------------------------//
/*!
 * Get all the state for the current level
 */
auto StateView::level_state() const -> const LevelState&
{
    CELER_EXPECT(state_);
    return state_.level_state[state_.level];
}

//---------------------------------------------------------------------------//
/*!
 * Get the current universe
 */
UniverseId StateView::universe() const
{
    const auto& cur_level = this->level_state();
    CELER_ENSURE(cur_level.univ_id);
    return cur_level.univ_id;
}

//---------------------------------------------------------------------------//
/*!
 * Get the current cell
 */
VolumeId StateView::volume() const
{
    const auto& cur_level = this->level_state();
    CELER_ENSURE(cur_level.cell);
    return cur_level.cell;
}

//---------------------------------------------------------------------------//
/*!
 * Get the current surface and universe it belongs to
 */
auto StateView::surface() const -> UnivSurface
{
    CELER_EXPECT(state_);
    UnivSurface result;

    if (state_.on_surface && state_.movement <= 0)
    {
        // On surface and hasn't moved
        CELER_ASSERT(state_.on_surface.level >= 0
                     && state_.on_surface.level <= state_.level);
        result.univ_id = state_.level_state[state_.on_surface.level].univ_id;
        result.surface = state_.on_surface.surface;
    }
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Whether a straight-line intercept has been found.
 */
bool StateView::has_next_distance() const
{
    CELER_ENSURE(state_.next_level <= state_.level);
    return state_.next_level >= 0;
}

//---------------------------------------------------------------------------//
/*!
 * Get the smallest straight-line distance to a boundary.
 *
 * This also accounts for any internal within-cell movement since the last
 * surface crossing or direction change.
 */
real_type StateView::next_distance() const
{
    CELER_EXPECT(this->has_next_distance());
    const Intersection& next = state_.level_state[state_.next_level].next;
    real_type           dist = next.distance - state_.movement;
    CELER_ENSURE(dist >= 0);
    return dist;
}

//---------------------------------------------------------------------------//
/*!
 * Get the next surface along the current direction.
 */
auto StateView::next_surface() const -> UnivSurface
{
    CELER_EXPECT(state_);
    UnivSurface result;

    if (this->has_next_distance())
    {
        CELER_ASSERT(state_.next_level <= state_.level);
        const LevelState& local = state_.level_state[state_.next_level];
        result.univ_id          = local.univ_id;
        result.surface          = local.next.surface;
    }
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Only the global universe has "move within boundary".
 *
 * \pre `intersect` or `distance_to_boundary` has been called.
 */
void StateView::move_internal(real_type distance)
{
    CELER_EXPECT(this->has_next_distance());
    CELER_EXPECT(distance > 0);
    CELER_EXPECT(distance <= this->next_distance());
    state_.movement += distance;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the current global position, accounting for local movement
 */
Real3 StateView::calc_position() const
{
    CELER_EXPECT(state_);

    // Result = global position + correction + movement * global direction
    Real3 result;
    for (int ax : range(3))
    {
        result[ax] = state_.level_state[0].pos[ax] + state_.correction[ax]
                     + state_.movement * state_.level_state[0].dir[ax];
    }
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Update position at all levels.
 *
 * This propagates the state's along-straight-line movement into the actual
 * multi-level positions (and distance-to-boundary), then clears the surface
 * flags (which may or may not be set).
 */
void StateView::apply_movement()
{
    CELER_EXPECT(state_.movement > 0);

    this->apply_movement_impl(state_.movement);

    // Reset movement
    state_.movement = 0;
    // No longer on a surface
    state_.on_surface = {};
}

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Update position at all levels.
 *
 * This propagates the state's along-straight-line movement into the actual
 * multi-level positions (and distance-to-boundary), then clears the surface
 * flags (which may or may not be set).
 */
void StateView::apply_movement_impl(real_type movement)
{
    CELER_EXPECT(state_);
    CELER_EXPECT(movement != 0);

    using axpy;
    using make_fixed_view;

    for (int level : range(state_.level + 1))
    {
        auto& ls = state_.level_state[level];
        CELER_ASSERT(ls.next.distance >= movement);
        ls.next.distance -= movement;
        axpy(movement, make_fixed_view(ls.dir), make_fixed_view(ls.pos));
    }
}

//---------------------------------------------------------------------------//
} // namespace celeritas
