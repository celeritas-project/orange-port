//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/LocalStateBuilder.hh
 * \brief LocalStateBuilder class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Assert.hh"
#include "../Definitions.hh"
#include "../TrackingState.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Build a result state reference.
 */
class LocalStateBuilder
{
  public:
    // Constructor
    explicit inline LocalStateBuilder(const TrackingState* state);

    // Get the state reference
    inline LocalState operator()(int level);

  private:
    //// DATA ////

    const TrackingState* state_;
};

//---------------------------------------------------------------------------//
// INLINE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Get the state reference.
 */
LocalStateBuilder::LocalStateBuilder(const TrackingState* state)
    : state_(state)
{
    CELER_EXPECT(state_);
}

//---------------------------------------------------------------------------//
/*!
 * Get the state reference.
 */
LocalState LocalStateBuilder::operator()(int level)
{
    CELER_EXPECT(level >= 0 && level <= state_->level);

    LocalState result;
    const auto&        level_state = state_->level_state[level];

    result.pos  = make_fixed_view(level_state.pos);
    result.dir  = make_fixed_view(level_state.dir);
    result.cell = level_state.cell;

    if (level == state_->on_surface.level)
    {
        // Deepest level uses surface state; others do not
        result.surface = state_->on_surface.surface;
        result.sense   = state_->on_surface.sense;
    }
    else
    {
        result.surface = {};
        result.sense   = {};
    }

    result.temp_senses    = &state_->temp_senses;
    result.temp_face_dist = &state_->temp_face_dist;
    return result;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
