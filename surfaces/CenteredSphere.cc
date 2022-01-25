
//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CenteredSphere.cc
 * \brief CenteredSphere class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CenteredSphere.hh"
#include <cmath>
#include <iostream>
#include "base/Assert.hh"
#include "orange/Transform.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "Sphere.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * Whether the given quadric can be simplified to this shape
 */
bool CenteredSphere::can_simplify(const Sphere& sq)
{
    const real_type tol = fuzziness().surface_simplification_abs();

    return (dot_product(sq.origin(), sq.origin()) <= ipow<2>(tol));
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * User-accessible construction
 */
CenteredSphere::CenteredSphere(real_type radius) : radius_sq_ipow<2>(radius)
{
    CELER_EXPECT(radius >= 0);
    CELER_ENSURE(radius_sq_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from a sphere
 */
CenteredSphere::CenteredSphere(const Sphere& sph) : radius_sq_(sph.radius_sq())
{
    CELER_EXPECT(vector_magnitude(sph.origin())
                 <= fuzziness().surface_simplification_abs());
    CELER_ENSURE(radius_sq_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * Return a translated sphere
 */
Sphere CenteredSphere::translated(const Transform& t) const
{
    return this->transformed(t);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a transformed sphere
 */
Sphere CenteredSphere::transformed(const Transform& t) const
{
    // Translation vector is the new origin of the sphere, ideally we could
    // just pass rsq instead of having to calculate the square root
    return Sphere(t.translation(), std::sqrt(radius_sq_));
}

//---------------------------------------------------------------------------//
/*!
 * Clip a bounding box to this surface
 */
void CenteredSphere::clip(Sense sense, BoundingBox& bbox) const
{
    if (sense == neg)
    {
        real_type radius = std::sqrt(radius_sq_);
        // Increase very slightly (bbox_rel is BVH fuzziness ~ 1e-7; we only
        // need to avoid clipping the actual sphere edges..)
        radius *= 1 + fuzziness().bbox_rel() * fuzziness().bbox_rel();
        for (int ax = 0; ax < 3; ++ax)
        {
            bbox.clip_lower(ax, -radius);
            bbox.clip_upper(ax, radius);
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * Output a very short description of this sufrace
 */
std::ostream& operator<<(std::ostream& os, const CenteredSphere& s)
{
    os << "Sphere: r=" << std::sqrt(s.radius_sq());
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
