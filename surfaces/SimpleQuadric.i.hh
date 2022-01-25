//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/SimpleQuadric.i.hh
 * \brief SimpleQuadric inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "base/Assert.hh"
#include "base/SoftEquivalence.hh"
#include "base/VectorFunctions.hh"
#include "detail/SolveQuadratic.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Inline construction from flattened coefficients (copying data).
 */
SimpleQuadric::SimpleQuadric(const real_type* coeff)
    : a_(coeff[0])
    , b_(coeff[1])
    , c_(coeff[2])
    , d_(coeff[3])
    , e_(coeff[4])
    , f_(coeff[5])
    , g_(coeff[6])
{
}

//---------------------------------------------------------------------------//
/*!
 * Determine the sense of the position relative to this surface (init).
 */
SignedSense SimpleQuadric::calc_sense(const Real3& pos) const
{
    const real_type x = pos[0];
    const real_type y = pos[1];
    const real_type z = pos[2];

    real_type result = (a_ * ipow<2>(x) + b_ * ipow<2>(y) + c_ * ipow<2>(z))
                       + (d_ * x + e_ * y + f_ * z) + (g_);

    return real_to_sense(result);
}

//---------------------------------------------------------------------------//
/*!
 * Determine intersection (numeric_limits::max() if no intersection).
 */
void SimpleQuadric::calc_intersections(const Real3& pos,
                                       const Real3& dir,
                                       bool         on_surface,
                                       real_type*   dist_iter) const
{
    CELER_EXPECT(dist_iter != nullptr);
    CELER_EXPECT(soft_equiv(1.0, vector_magnitude(dir), 1.0e-6));

    real_type* const end_iter = dist_iter + num_intersections();

    // Translated position
    const real_type x = pos[0];
    const real_type u = dir[0];
    const real_type y = pos[1];
    const real_type v = dir[1];
    const real_type z = pos[2];
    const real_type w = dir[2];

    // Quadratic values
    real_type a      = (a_ * u) * u + (b_ * v) * v + (c_ * w) * w;
    real_type half_b = (a_ * x + half * d_) * u + (b_ * y + half * e_) * v
                       + (c_ * z + half * f_) * w;
    real_type c = (a_ * x + d_) * x + (b_ * y + e_) * y + (c_ * z + f_) * z
                  + g_;

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
Real3 SimpleQuadric::calc_normal(const Real3& pos) const
{
    const real_type x = pos[0];
    const real_type y = pos[1];
    const real_type z = pos[2];

    Real3 norm;
    norm[0] = a_ * x + d_;
    norm[1] = b_ * y + e_;
    norm[2] = c_ * z + f_;

    normalize_direction(&norm);
    return norm;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
