//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CenteredSphere.i.hh
 * \brief CenteredSphere inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "base/VectorFunctions.hh"
#include "detail/SolveQuadratic.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct from coeffients
 */
CenteredSphere::CenteredSphere(const real_type* coeff) : radius_sq_(coeff[0])
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine the sense of the position relative to this surface (init)
 */
SignedSense CenteredSphere::calc_sense(const Real3& pos) const
{
    using Axis::x;
    using Axis::y;
    using Axis::z;
    const real_type x = pos[X];
    const real_type y = pos[Y];
    const real_type z = pos[Z];

    return real_to_sense(ipow<2>(x) + ipow<2>(y) + ipow<2>(z) - radius_sq_);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine the distance to intersection with this sphere
 */
void CenteredSphere::calc_intersections(const Real3& pos,
                                        const Real3& dir,
                                        bool         on_surface,
                                        real_type*   dist_iter) const
{
    CELER_EXPECT(dist_iter != nullptr);
    using Axis::x;
    using Axis::y;
    using Axis::z;

    real_type* const end_iter = dist_iter + num_intersections();

    const real_type x = pos[X];
    const real_type y = pos[Y];
    const real_type z = pos[Z];

    // (X - X0) .. U
    real_type half_b = x * dir[X] + y * dir[Y] + z * dir[Z];

    if (!on_surface)
    {
        // (X - X0) . (X - X0) - R^2
        real_type c = ipow<2>(x) + ipow<2>(y) + ipow<2>(z) - radius_sq_;
        dist_iter   = detail::solve_quadratic(half_b, c, dist_iter);
    }
    else
    {
        dist_iter = detail::solve_degenerate_quadratic1(half_b, dist_iter);
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
Real3 CenteredSphere::calc_normal(const Real3& pos) const
{
    Real3 norm = pos.vector();

    normalize_direction(&norm);
    return norm;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
