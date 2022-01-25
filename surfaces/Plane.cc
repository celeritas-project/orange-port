//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/Plane.cc
 * \brief Plane class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Plane.hh"

#include <cmath>
#include <iostream>
#include "base/FixedViewArray.hh"
#include "base/VectorFunctions.hh"
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
bool Plane::can_simplify(const SimpleQuadric& sq)
{
    auto            second = sq.second();
    const real_type tol    = fuzziness().surface_simplification_abs();
    return std::fabs(second[0]) < tol && std::fabs(second[1]) < tol
           && std::fabs(second[2]) < tol;
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * Construct with unnormalized normal and point on the plane.
 */
Plane::Plane(const Real3& n, const Real3& p) : normal_(make_vector(n))
{
    CELER_EXPECT(std::fabs(vector_magnitude(normal_)) > 0.0);

    // Normalize the vector to make it a unit-normal
    normalize_direction(&normal_);

    // Displacement is the dot product of the point and the normal
    d_ = dot_product(normal_, p);
}

//---------------------------------------------------------------------------//
/*!
 * Construct with normalized normal and displacement.
 *
 * The displacement along the normal vector will be normalized according to the
 * normal.
 */
Plane::Plane(const Real3& n, real_type d) : normal_(make_vector(n)), d_(d)
{
    CELER_EXPECT(soft_equiv(vector_magnitude(n),
                            1.0,
                            100 * fuzziness().surface_simplification_abs()));

    // Calculate magnitude of the normal if it's not quite 1
    auto norm = vector_magnitude(normal_);
    CELER_ASSERT(norm > 0);
    normal_ /= norm;

    // Renormalize intersection point based on the new normal (preserve
    // intersection, not displacement)
    d_ /= norm;
}

//---------------------------------------------------------------------------//
/*!
 * Construct from a degenerate quadric
 */
Plane::Plane(const SimpleQuadric& sq)
{
    CELER_EXPECT(can_simplify(sq));

    // Normal is same as first-order coefficients
    normal_ = make_vector(sq.first());

    // Calculate normalization for normal component, since the SQ doesn't have
    // to be scaled
    auto norm = vector_magnitude(normal_);
    CELER_ASSERT(norm > 0);
    normal_ /= norm;

    // Simple quadric is calculated as dx + ey + fz + g  = 0, but
    // plane is ax + bx + cz - d = 0, so we need to reverse the sign of the
    // scalar component. We also need to divide by the normal component.
    d_ = -sq.zeroth() / norm;

    CELER_ENSURE(soft_equiv(vector_magnitude(normal_), 1.0));
}

//---------------------------------------------------------------------------//
/*!
 * Return a new plane translated by some vector
 */
Plane Plane::translated(const Transform& t) const
{
    CELER_EXPECT(!t.has_rotation());

    Plane result(*this);

    real_type displacement = dot_product(t.translation(), normal());
    result.d_ += displacement;
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Return a plane transformed by a rotate/translate
 */
Plane Plane::transformed(const Transform& t) const
{
    CELER_EXPECT(t.has_rotation());

    // Calculate new normal
    Real3 rn(normal());
    t.rotate_to_parent(rn);

    // Get the original point on the plane
    Real3 rp(normal());
    rp *= d_;

    // Transform the point
    t.daughter_to_parent(rp);

    // Build the general plane corresponding to this transform
    return Plane(rn, dot_product(rn, rp));
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Describe to a stream.
 */
std::ostream& operator<<(std::ostream& os, const Plane& s)
{
    os << "Plane: n=(" << s.normal() << "), d=" << s.displacement();
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
