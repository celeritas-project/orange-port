//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/Sphere.cc
 * \brief Sphere class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Sphere.hh"

#include <cmath>
#include <iostream>
#include "base/FixedViewArray.hh"
#include "base/Assert.hh"
#include "orange/Transform.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "SimpleQuadric.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * Whether the given quadric can be simplified to this shape
 */
bool Sphere::can_simplify(const SimpleQuadric& sq)
{
    // TODO: implement this
    return false;
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief User-accessible construction
 */
Sphere::Sphere(const Real3& origin, real_type radius)
    : origin_(make_vector(origin)), radius_sq_ipow<2>(radius)
{
    CELER_EXPECT(radius >= 0);

    CELER_ENSURE(radius_sq_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from a degenerate SQ
 */
Sphere::Sphere(const SimpleQuadric& sq)
{
    NotImplemented("simplifying a simple quadric to a sphere");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a new sphere transformed
 */
Sphere Sphere::transformed(const Transform& t) const
{
    // Translate/rotate sphere's origin; rotation of sphere itself is obviously
    // invariant
    Sphere s(*this);
    t.daughter_to_parent(s.origin_);
    return s;
}

//---------------------------------------------------------------------------//
/*!
 * Clip a bounding box to this surface
 */
void Sphere::clip(Sense sense, BoundingBox& bbox) const
{
    if (sense == neg)
    {
        real_type radius = std::sqrt(radius_sq_);
        // Increase very slightly (bbox_rel is BVH fuzziness ~ 1e-7; we only
        // need to avoid clipping the actual sphere edges..)
        radius *= 1 + fuzziness().bbox_rel() * fuzziness().bbox_rel();
        for (int ax = 0; ax < 3; ++ax)
        {
            bbox.clip_lower(ax, origin_[ax] - radius);
            bbox.clip_upper(ax, origin_[ax] + radius);
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * Output a very short description of this sufrace
 */
std::ostream& operator<<(std::ostream& os, const Sphere& s)
{
    os << "Sphere: r=" << std::sqrt(s.radius_sq()) << " at " << s.origin();
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
