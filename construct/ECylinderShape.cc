//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ECylinderShape.cc
 * \brief ECylinderShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ECylinderShape.hh"

#include <algorithm>
#include <cmath>

#include "base/Assert.hh"
#include "base/Constants.hh"
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

using Axis::x;
using Axis::y;

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Constructor.
 */
ECylinderShape::ECylinderShape(PlaneVector radii, real_type lo_z, real_type hi_z)
    : radii_(radii), lo_(lo_z), hi_(hi_z)
{
    for (Axis ax : {X, Y})
    {
        CELER_VALIDATE(radii_[ax] >= 0,
                       << "Elliptical cylinder radius along XY"[ax]
                       << " must be positive (given r=" << radii_[ax] << ")");
    }
    CELER_VALIDATE(lo_ <= hi_,
                   << "Lower extent " << lo_
                   << " must be less than upper extent " << hi_
                   << " for elliptical cylinder");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string ECylinderShape::type() const
{
    return "ecylinder";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool ECylinderShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief ECylinderShape volume
 */
real_type ECylinderShape::volume() const
{
    using constants::pi;

    real_type vol = pi * radii_[X] * radii_[Y] * (hi_ - lo_);
    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type ECylinderShape::inradius() const
{
    real_type max_radius = std::max(radii_[X], radii_[Y]);
    real_type hheight    = half * (hi_ - lo_);
    return std::sqrt(ipow<2>(hheight) + 2.0 * ipow<2>(max_radius));
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void ECylinderShape::build(ShapeBuilder& build) const
{
    // Second-order coefficients are product of the other squared radius;
    // Zeroth-order coefficient is the product of both squared radii
    PlaneVector rsq = ipow<2>(radii_);

    const Real3     abc(rsq[Y], rsq[X], 0);
    const real_type g = -rsq[X] * rsq[Y];
    build.simple_quadric(abc, Real3(0, 0, 0), g, Real3(0, 0, 0));

    // Build the end caps
    build.plane(Axis::z, pos, lo_);
    build.plane(Axis::z, neg, hi_);

    // Clip radial directions based on outer shape radius
    for (Axis ax : {X, Y})
    {
        build.clip_bbox_lower(ax, -radii_[ax]);
        build.clip_bbox_upper(ax, radii_[ax]);
    }
}

//---------------------------------------------------------------------------//
} // namespace celeritas
