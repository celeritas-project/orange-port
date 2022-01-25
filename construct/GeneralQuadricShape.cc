//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/GeneralQuadricShape.cc
 * \brief GeneralQuadricShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "GeneralQuadricShape.hh"

#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
GeneralQuadricShape::GeneralQuadricShape(real_type a,
                                         real_type b,
                                         real_type c,
                                         real_type d,
                                         real_type e,
                                         real_type f,
                                         real_type g,
                                         real_type h,
                                         real_type i,
                                         real_type j)
    : a_(a), b_(b), c_(c), d_(d), e_(e), f_(f), g_(g), h_(h), i_(i), j_(j)
{
    CELER_VALIDATE(
        std::any_of(&a_, &j_ + 1, [](real_type v) { return v != 0; }),
        << "At least one quadric coefficient must be nonzero");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string GeneralQuadricShape::type() const
{
    return "general_quadric";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings).
 *
 * Cones and other hyperbolic structures have two separated components and thus
 * are not necessarily convex. Theoretically we could figure out whether the
 * quadric here is a two-parter or just a single continuous convex surface
 * (e.g. ellipsoid), perhaps by using the Surface simplification logic, but
 * it's probably not worth the effort, so we just assume the worst.
 */
bool GeneralQuadricShape::is_convex() const
{
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * Shape volume
 */
real_type GeneralQuadricShape::volume() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type GeneralQuadricShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void GeneralQuadricShape::build(ShapeBuilder& build) const
{
    build.general_quadric({a_, b_, c_}, {d_, e_, f_}, {g_, h_, i_}, j_, neg);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
