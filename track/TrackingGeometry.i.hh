//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/TrackingGeometry.i.hh
 * \brief TrackingGeometry inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "StateView.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// ACCESSORS
//---------------------------------------------------------------------------//
/*!
 * Get the universe for this ID
 */
const Universe& TrackingGeometry::get(UniverseId univ_id) const
{
    CELER_EXPECT(univ_id < universes_.size());
    return *universes_[univ_id.get()];
}

//---------------------------------------------------------------------------//
// TRACKING
//---------------------------------------------------------------------------//
/*!
 * Get or calculate the distance-to-boundary
 */
real_type TrackingGeometry::distance_to_boundary(State_t& state) const
{
    CELER_EXPECT(state);
    StateView view(state);

    if (!view.has_next_distance())
    {
        // Calculate the distance-to-boundary
        CELER_ASSERT(state.movement <= 0);
        this->intersect(state);
        CELER_ASSERT(view.has_next_distance());
    }

    // Calculated next-distance minus any straight-line movement
    real_type distance = view.next_distance();
    CELER_ENSURE(distance >= 0);
    return distance;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
