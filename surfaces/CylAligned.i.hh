//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CylAligned.i.hh
 * \brief CylAligned inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "base/Assert.hh"
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
constexpr auto CylAligned<T>::surface_type() -> SurfaceType
{
    return (T == Axis::X ? SurfaceType::cx
                         : (T == Axis::Y ? SurfaceType::cy : SurfaceType::cz));
}

//---------------------------------------------------------------------------//
// INLINE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Construct from coeffients
 */
template<Axis T>
CylAligned<T>::CylAligned(const real_type* coeff)
    : origin_u_(coeff[0]), origin_v_(coeff[1]), radius_sq_(coeff[2])
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine the sense of the position relative to this surface (init)
 */
template<Axis T>
SignedSense CylAligned<T>::calc_sense(const Real3& pos) const
{
    const real_type u = pos[u_index()] - origin_u_;
    const real_type v = pos[v_index()] - origin_v_;

    return real_to_sense(ipow<2>(u) + ipow<2>(v) - radius_sq_);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine distance to intersect.
 */
template<Axis T>
void CylAligned<T>::calc_intersections(const Real3& pos,
                                       const Real3& dir,
                                       bool         on_surface,
                                       real_type*   dist_iter) const
{
    CELER_EXPECT(dist_iter != nullptr);

    real_type* const end_iter = dist_iter + num_intersections();

    // 1 - \omega \dot e
    const real_type a = 1 - dir[t_index()] * dir[t_index()];

    // No intersection if we're traveling along the cylinder axis
    if (a)
    {
        const real_type u = pos[u_index()] - origin_u_;
        const real_type v = pos[v_index()] - origin_v_;

        const real_type a_inv = 1. / a;

        // \omega \dot (x - x_0)
        real_type half_b = dir[u_index()] * u + dir[v_index()] * v;
        half_b *= a_inv;

        if (!on_surface)
        {
            // (x - x_0) \dot (x - x_0) - ipow<2>(R)
            real_type c = a_inv * (ipow<2>(u) + ipow<2>(v) - radius_sq_);
            dist_iter   = detail::solve_quadratic(half_b, c, dist_iter);
        }
        else
        {
            dist_iter = detail::solve_degenerate_quadratic1(half_b, dist_iter);
        }
    }

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
Real3 CylAligned<T>::calc_normal(const Real3& pos) const
{
    Real3 norm(0);

    norm[u_index()] = pos[u_index()] - origin_u_;
    norm[v_index()] = pos[v_index()] - origin_v_;

    normalize_direction(&norm);
    return norm;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
