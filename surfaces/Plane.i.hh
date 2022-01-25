//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/Plane.i.hh
 * \brief Plane inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "base/Assert.hh"
#include "base/SoftEquivalence.hh"
#include "base/VectorFunctions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Inline construction from flattened coefficients (copying data).
 */
Plane::Plane(const real_type* coeff)
    : normal_(coeff[0], coeff[1], coeff[2]), d_(coeff[3])
{
}

//---------------------------------------------------------------------------//
/*!
 * Determine the sense of the position relative to this surface (init).
 */
SignedSense Plane::calc_sense(const Real3& x) const
{
    return real_to_sense(dot_product(normal_, x) - d_);
}

//---------------------------------------------------------------------------//
/*!
 * Determine intersection (numeric_limits::max() if no intersection).
 */
void Plane::calc_intersections(const Real3& pos,
                               const Real3& dir,
                               bool         on_surface,
                               real_type*   dist_iter) const
{
    CELER_EXPECT(dist_iter != nullptr);
    CELER_EXPECT(soft_equiv(1.0, vector_magnitude(dir), 1.0e-6));

    // Calculate n\dot\Omega
    real_type n_dot_omega = dot_product(normal_, dir);

    if (!on_surface && n_dot_omega != 0)
    {
        *dist_iter = (d_ - dot_product(normal_, pos)) / n_dot_omega;
        if (*dist_iter < 0)
        {
            // Past the surface
            *dist_iter = no_intersection();
        }
    }
    else
    {
        // It's on or parallel to the surface
        *dist_iter = no_intersection();
    }
}

//---------------------------------------------------------------------------//
} // namespace celeritas
