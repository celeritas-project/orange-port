//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/StateView.hh
 * \brief StateView class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <utility>
#include "TrackingState.hh"
#include "../Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Utility functions for accessing the tracking state.
 *
 * This provides a level of abstraction for accessing state data. It's modeled
 * after the Store/View paradigm in Celeritas to provide a separation between
 * data managment and computation, and it thus could be a first step toward a
 * GPU implementation.
 */
class StateView
{
  public:
    //@{
    //! Public type aliases
    using LevelState = TrackingState::LevelState;
    //@}

    //! Basic information needed to initialize the state
    struct Initializer_t
    {
        SpanConstReal3 pos;
        SpanConstReal3 dir;
    };

    //! Whether and where the state is on a surface
    struct UnivSurface
    {
        UniverseId univ_id;
        SurfaceId  surface;

        operator bool() const { return bool(univ_id); }
    };

  public:
    // Constructor
    inline explicit StateView(const TrackingState& state);

    // Initialize the state
    inline StateView& operator=(const Initializer_t&);

    //// ACCESSORS ////

    // State at the highest level
    inline const LevelState& global_state() const;

    // State at the current level
    inline const LevelState& level_state() const;

    // Get the current universe
    inline UniverseId universe() const;

    // Get the current cell
    inline VolumeId volume() const;

    // Get the current surface (and universe ID), may be false
    inline UnivSurface surface() const;

    // Whether a straight-line intercept has been found
    inline bool has_next_distance() const;

    // Distance to next boundary along a straight line
    inline real_type next_distance() const;

    // Get the next surface (and universe ID), may be false
    inline UnivSurface next_surface() const;

    //// MOVEMENT ////

    // Only the global universe has "move within boundary"
    inline void move_internal(real_type distance);

    // Propagate 'movement' into actual multi-level positions
    inline void apply_movement();

    // Calculate the current global position, accounting for local movement
    inline Real3 calc_position() const;

  private:
    //// DATA ////

    // Mutable, but not changed by 'const' methods
    TrackingState& state_;

    //// HELPER FUNCTIONS ////

    // Propagate 'movement' into actual multi-level positions
    inline void apply_movement_impl(real_type movement);
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "StateView.i.hh"
//---------------------------------------------------------------------------//
