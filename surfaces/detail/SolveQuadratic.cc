//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/detail/SolveQuadratic.cc
 * \brief Quadratic solvers
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include "base/Assert.hh"
#include "orange/Fuzziness.hh"
#include "SolveQuadratic.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Find all nonnegative roots for general quadric surfaces.
 *
 * This is used for cones, simple quadrics, and general quadrics.
 */
real_type* solve_quadratic_general(SurfaceState on_surface,
                                   real_type    a,
                                   real_type    half_b,
                                   real_type    c,
                                   real_type*   x)
{
    // Calculate 1/a, accounting for along_surface values
    bool      along_surface = std::fabs(a) < Fuzziness::quadratic_abs();
    real_type a_inv = 1.0 / (!along_surface ? a : Fuzziness::quadratic_abs());

    // Rescale b/2 and c so that a = 1 in the quadratic equation
    half_b *= a_inv;
    c *= a_inv;

    if (!on_surface && !along_surface)
    {
        return solve_quadratic(half_b, c, x);
    }
    else if (!on_surface && along_surface)
    {
        return solve_degenerate_quadratic2(half_b, c, x);
    }
    else if (on_surface && !along_surface)
    {
        return solve_degenerate_quadratic1(half_b, x);
    }

    // On surface *and* along it: no intersection
    return x;
}

//---------------------------------------------------------------------------//
/*!
 * Find all nonnegative roots of x^2 + b*x + c = 0.
 *
 * Callers:
 * - General quadratic solve: not on nor along surface
 * - Sphere when not on surface
 * - Cylinder when not on surface
 *
 * \param[in]  half_b  b/2
 * \param[in]  c    c
 * \param[out] x    roots, if any, in order of increasing magnitude
 *
 * \return Pointer to past-the-end root
 */
real_type* solve_quadratic(real_type half_b, real_type c, real_type* x)
{
    // In this case, the quadratic formula can be written as:
    // x = -b/2 +/- [(b/2)^2 - c]^half.

    real_type b2_4 = ipow<2>(half_b); // (b/2)^2

    if (b2_4 < c)
    {
        // No real roots
    }
    else if (b2_4 > c)
    {
        // Two real roots, r1 and r2
        real_type t2 = std::sqrt(b2_4 - c); // (b^2 - 4ac) / 4
        real_type r1 = -half_b - t2;
        real_type r2 = -half_b + t2; // r2 > r1

        if (r1 >= 0)
        {
            *x++ = r1;
            *x++ = r2;
        }
        else if (r2 >= 0)
        {
            *x++ = r2;
        }
    }
    else
    {
        // One real root, r1
        real_type r1 = -half_b;

        if (r1 >= 0)
        {
            *x++ = r1;
        }
    }

    return x;
}

//---------------------------------------------------------------------------//
/*!
 * Find positive root of x^2 + b*x + c = 0, where x = 0 is a root.
 *
 * Callers:
 * - General quadratic solve above (on but not along surface)
 * - Sphere when on surface
 * - Cylinder when on surface
 *
 * \param[in]  half_b  b/2
 * \param[out] x    root, if any
 *
 * \return Pointer to past-the-end root
 */
real_type* solve_degenerate_quadratic1(real_type half_b, real_type* x)
{
    // If x = 0 is a root, then c = 0 and x = -b is the other root
    real_type r2 = -2.0 * half_b;

    if (r2 > 0)
    {
        *x++ = r2;
    }

    return x;
}

//---------------------------------------------------------------------------//
/*!
 * Find nonnegative roots of a*x^2 + b*x + c = 0, where a -> 0.
 *
 * Callers:
 * - General quadratic solve above (along but not on the surface)
 *
 * \param[in]  b_a_2  (b/a)/2
 * \param[in]  c_a    c/a
 * \param[out] x      positive root, if any
 *
 * \return Pointer to past-the-end root
 *
 * Note that as both a and b -> 0, |x| -> infinity. Thus, if b/a is
 * sufficiently small, no positive root is returned (because this case
 * corresponds to a ray crossing a surface at an extreme distance).
 */
real_type*
solve_degenerate_quadratic2(real_type b_a_2, real_type c_a, real_type* x)
{
    // If a -> 0, then x = -c/b
    if (std::fabs(b_a_2) > Fuzziness::quadratic_abs())
    {
        real_type r1 = -c_a / (2.0 * b_a_2);

        if (r1 >= 0)
        {
            *x++ = r1;
        }
    }

    // As a and b -> 0, |x| -> infinity
    return x;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
