//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CylinderShape.cc
 * \brief CylinderShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CylinderShape.hh"

#include <cmath>

#include "base/Assert.hh"
#include "base/Constants.hh"
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct the cylinder
 */
CylinderShape::CylinderShape(Axis      axis,
                             real_type radius,
                             real_type lo,
                             real_type hi)
    : axis_(axis), radius_(radius), lo_(lo), hi_(hi)
{
    CELER_VALIDATE(axis < 3, << "Invalid axis " << axis << " for cylinder");
    CELER_VALIDATE(lo <= hi,
                   << "Lower extent " << lo
                   << " must be less than upper extent " << hi
                   << " for cylinder");
    CELER_VALIDATE(radius > 0, << "Radius for cylinder must be positive");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct the cylinder
 */
CylinderShape::CylinderShape(real_type radius, real_type height)
    : axis_(Axis::z), radius_(radius), lo_(-height / 2), hi_(height / 2)
{
    CELER_VALIDATE(height > 0, << "Negative cylinder height " << height);
    CELER_VALIDATE(radius > 0, << "Negative cylinder radius " << radius);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Shape class name
 */
std::string CylinderShape::type() const
{
    return "cylinder";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool CylinderShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief CylinderShape volume
 */
real_type CylinderShape::volume() const
{
    using constants::pi;

    return pi * ipow<2>(radius_) * (hi_ - lo_);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Largest sphere radius that fits in this shape
 */
real_type CylinderShape::inradius() const
{
    return std::min(radius_, (hi_ - lo_) / 2);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct surfaces for this shape
 */
void CylinderShape::build(ShapeBuilder& build) const
{
    // Build cylinder
    build.cyl(axis_, neg, radius_);

    // Build bounding planes (above lower, below upper)
    build.plane(axis_, pos, lo_);
    build.plane(axis_, neg, hi_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
