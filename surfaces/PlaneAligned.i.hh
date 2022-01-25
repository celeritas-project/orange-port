//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/PlaneAligned.i.hh
 * \brief PlaneAligned method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "base/Assert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// INLINE STATIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Get the surface type of this surface.
 */
template<Axis T>
constexpr auto PlaneAligned<T>::surface_type() -> SurfaceType
{
    return (T == Axis::X ? SurfaceType::px
                         : (T == Axis::Y ? SurfaceType::py : SurfaceType::pz));
}

//---------------------------------------------------------------------------//
// INLINE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Construct from coeffients
 */
template<Axis T>
PlaneAligned<T>::PlaneAligned(const real_type* coeff) : position_(coeff[0])
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine the sense of the position relative to this surface (init)
 */
template<Axis T>
SignedSense PlaneAligned<T>::calc_sense(const Real3& x) const
{
    return real_to_sense(x[t_index()] - position_);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine distance to intersect.
 *
 * This will set the value to "no_intersection()" if it never intersects.
 */
template<Axis T>
void PlaneAligned<T>::calc_intersections(const Real3& pos,
                                         const Real3& dir,
                                         bool         on_surface,
                                         real_type*   dist_iter) const
{
    CELER_EXPECT(dist_iter != nullptr);

    if (!on_surface && dir[t_index()] != 0)
    {
        *dist_iter = (position_ - pos[t_index()]) / dir[t_index()];
        if (*dist_iter < 0)
        {
            // Past the surface
            *dist_iter = no_intersection();
        }
    }
    else
    {
        // If it's heading away from, on, or parallel to the surface
        *dist_iter = no_intersection();
    }
}

//---------------------------------------------------------------------------//
/*!
 * Calculate outward normal at a position
 */
template<Axis T>
Real3 PlaneAligned<T>::calc_normal(const Real3&) const
{
    Real3 norm(0);
    norm[t_index()] = 1..;
    return norm;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
