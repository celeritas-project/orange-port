//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/SlabShape.cc
 * \brief SlabShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SlabShape.hh"

#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct the SlabShape
 */
SlabShape::SlabShape(Axis axis, real_type lo, real_type hi)
    : axis_(axis), lo_(lo), hi_(hi)
{
    CELER_VALIDATE(axis < 3, << "Invalid axis " << axis << " for slab");
    CELER_VALIDATE(lo <= hi,
                   << "Lower extent " << lo
                   << " must be less than upper extent " << hi << " for slab");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Shape class name
 */
std::string SlabShape::type() const
{
    return "slab";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool SlabShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief SlabShape volume
 */
real_type SlabShape::volume() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type SlabShape::inradius() const
{
    return half * (hi_ - lo_);
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void SlabShape::build(ShapeBuilder& build) const
{
    // Build surfaces
    build.plane(axis_, pos, lo_);
    build.plane(axis_, neg, hi_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
