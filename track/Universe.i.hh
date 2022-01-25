//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/Universe.i.hh
 * \brief Universe inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Find the starting location for a particle *or* cross a surface.
 */
Initialization Universe::initialize(LocalState local_state) const
{
    return tracker_->initialize(local_state);
}

//---------------------------------------------------------------------------//
/*!
 * Find the distance-to-intercept at this level.
 */
Intersection Universe::intersect(LocalState local_state) const
{
    return tracker_->intersect(local_state);
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the local surface normal.
 */
Real3 Universe::normal(LocalState local_state) const
{
    return tracker_->normal(local_state);
}

//---------------------------------------------------------------------------//
/*!
 * Get the local daughter for the given cell or nullptr if not present.
 */
auto Universe::daughter(VolumeId cell) const -> const Daughter*
{
    CELER_EXPECT(cell < this->num_volumes());

    if (daughters_.empty())
    {
        return nullptr;
    }

    const Daughter& d = daughters_[cell.unchecked_get()];
    if (!d.universe)
    {
        return nullptr;
    }
    return &d;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
