//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/ConeAligned.i.hh
 * \brief ConeAligned inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include <cmath>
#include "base/VectorFunctions.hh"
#include "detail/SolveQuadratic.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// INLINE STATIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Get the surface type of this surface.
 */
template<Axis T>
constexpr auto ConeAligned<T>::surface_type() -> SurfaceType
{
    return (T == Axis::X ? SurfaceType::kx
                         : (T == Axis::Y ? SurfaceType::ky : SurfaceType::kz));
}

//---------------------------------------------------------------------------//
// INLINE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Construct from coeffients
 */
template<Axis T>
ConeAligned<T>::ConeAligned(const real_type* coeff)
    : origin_(coeff[0], coeff[1], coeff[2]), tsq_(coeff[3])
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine the sense of the position relative to this surface (init)
 */
template<Axis T>
SignedSense ConeAligned<T>::calc_sense(const Real3& pos) const
{
    real_type x = pos[t_index()] - origin_[t_index()];
    real_type y = pos[u_index()] - origin_[u_index()];
    real_type z = pos[v_index()] - origin_[v_index()];

    return real_to_sense((-tsq_ * ipow<2>(x)) + ipow<2>(y) + ipow<2>(z));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine distance to intersect.
 *
 * \f[
    (y - yc)^2 + (z - zc)^2 - t^2 * (x - xc)^2 = 0
   \f]
 */
template<Axis T>
void ConeAligned<T>::calc_intersections(const Real3& pos,
                                        const Real3& dir,
                                        bool         on_surface,
                                        real_type*   dist_iter) const
{
    CELER_EXPECT(dist_iter != nullptr);

    real_type* const end_iter = dist_iter + num_intersections();

    // Expand translated positions into 'xyz' coordinate system
    const real_type x = pos[t_index()] - origin_[t_index()];
    const real_type y = pos[u_index()] - origin_[u_index()];
    const real_type z = pos[v_index()] - origin_[v_index()];

    const real_type u = dir[t_index()];
    const real_type v = dir[u_index()];
    const real_type w = dir[v_index()];

    // Scaled direction
    real_type a      = (-tsq_ * ipow<2>(u)) + ipow<2>(v) + ipow<2>(w);
    real_type half_b = (-tsq_ * x * u) + (y * v) + (z * w);
    real_type c      = (-tsq_ * ipow<2>(x)) + ipow<2>(y) + ipow<2>(z);

    // Solve quadric
    dist_iter
        = detail::solve_quadratic_general(on_surface, a, half_b, c, dist_iter);

    // Unrolled fill, maximum of 2 intersections
    if (dist_iter != end_iter)
        *dist_iter++ = no_intersection();
    if (dist_iter != end_iter)
        *dist_iter++ = no_intersection();
}

//---------------------------------------------------------------------------//
/*!
 * Calculate outward normal at a position
 */
template<Axis T>
Real3 ConeAligned<T>::calc_normal(const Real3& pos) const
{
    Real3 norm = pos.vector();
    norm -= origin_;
    norm[t_index()] *= -tsq_;

    normalize_direction(&norm);
    return norm;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
