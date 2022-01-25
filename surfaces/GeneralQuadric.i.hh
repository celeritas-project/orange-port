//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/GeneralQuadric.i.hh
 * \brief GeneralQuadric inline method definitions
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
GeneralQuadric::GeneralQuadric(const real_type* coeff)
    : a_(coeff[0])
    , b_(coeff[1])
    , c_(coeff[2])
    , d_(coeff[3])
    , e_(coeff[4])
    , f_(coeff[5])
    , g_(coeff[6])
    , h_(coeff[7])
    , i_(coeff[8])
    , j_(coeff[9])
{
}

//---------------------------------------------------------------------------//
/*!
 * Determine the sense of the position relative to this surface (init).
 */
SignedSense GeneralQuadric::calc_sense(const Real3& pos) const
{
    const real_type x = pos[0];
    const real_type y = pos[1];
    const real_type z = pos[2];

    real_type result = (a_ * x + d_ * y + f_ * z + g_) * x
                       + (b_ * y + e_ * z + h_) * y + (c_ * z + i_) * z + j_;

    return real_to_sense(result);
}

//---------------------------------------------------------------------------//
/*!
 * Determine intersection (numeric_limits::max() if no intersection).
 */
void GeneralQuadric::calc_intersections(const Real3& pos,
                                        const Real3& dir,
                                        bool         on_surface,
                                        real_type*   dist_iter) const
{
    CELER_EXPECT(dist_iter != nullptr);

    real_type* const end_iter = dist_iter + num_intersections();

    const real_type x = pos[0];
    const real_type u = dir[0];
    const real_type y = pos[1];
    const real_type v = dir[1];
    const real_type z = pos[2];
    const real_type w = dir[2];

    // Quadratic values
    real_type a = (a_ * u + d_ * v) * u + (b_ * v + e_ * w) * v
                  + (c_ * w + f_ * u) * w;
    real_type half_b = half
                       * ((2 * a_ * x + d_ * y + f_ * z + g_) * u
                          + (2 * b_ * y + d_ * x + e_ * z + h_) * v
                          + (2 * c_ * z + +e_ * y + f_ * x + i_) * w);
    real_type c = ((a_ * x + d_ * y + g_) * x + (b_ * y + e_ * z + h_) * y
                   + (c_ * z + f_ * x + i_) * z + j_);

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
Real3 GeneralQuadric::calc_normal(const Real3& pos) const
{
    const real_type x = pos[0];
    const real_type y = pos[1];
    const real_type z = pos[2];

    Real3 norm;
    norm[0] = 2.0 * a_ * x + d_ * y + f_ * z + g_;
    norm[1] = 2.0 * b_ * y + d_ * x + e_ * z + h_;
    norm[2] = 2.0 * c_ * z + e_ * y + f_ * x + i_;

    normalize_direction(&norm);
    return norm;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
