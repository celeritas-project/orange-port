//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ConeShape.cc
 * \brief ConeShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ConeShape.hh"

#include <algorithm>
#include <cmath>

#include "orange/Definitions.hh"
#include "base/Constants.hh"
#include "base/GeometryUtils.hh"
#include "base/Range.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Constructor.
 */
ConeShape::ConeShape(Axis axis, PairDbl radii, PairDbl extents)
    : axis_(axis)
    , lo_radius_(radii.first)
    , hi_radius_(radii.second)
    , lo_(extents.first)
    , hi_(extents.second)
{
    CELER_VALIDATE(axis_ < 3, << "Invalid axis " << axis << " for cone");
    CELER_VALIDATE(lo_radius_ >= 0,
                   << "Lower radius for cone must be nonnegative");
    CELER_VALIDATE(hi_radius_ >= 0,
                   << "Upper radius for cone must be nonnegative");
    CELER_VALIDATE(hi_radius_ != 0 || lo_radius_ != 0,
                   << "Can't build a cone shape with two zero radii");
    CELER_VALIDATE(lo_ <= hi_,
                   << "Lower extent " << lo_
                   << " must be less than upper extent " << hi_ << " for cone");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string ConeShape::type() const
{
    return "cone";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool ConeShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief ConeShape volume
 */
real_type ConeShape::volume() const
{
    real_type vol = constants::pi / 3 * (hi_ - lo_)
                    * (ipow<2>(lo_radius_) + lo_radius_ * hi_radius_
                       + ipow<2>(hi_radius_));
    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type ConeShape::inradius() const
{
    // Defined as the distance from the cone's origin
    // to the largest radii AABB vertex
    real_type hheight = half * (hi_ - lo_);
    real_type hwidth  = std::max(lo_radius_, hi_radius_);

    // Sqrt of the sum of the squares h^2 + w^2 + w^2
    return std::sqrt(ipow<2>(hheight) + 2.0 * ipow<2>(hwidth));
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void ConeShape::build(ShapeBuilder& build) const
{
    const real_type dr = std::fabs(lo_radius_ - hi_radius_);
    const real_type dz = hi_ - lo_;
    CELER_ASSERT(dz > 0);

    // Calculate vanishing point (origin)
    real_type origin_pos = 0;
    if (lo_radius_ > hi_radius_)
    {
        // Cone opens downward
        origin_pos = lo_ + lo_radius_ * dz / dr;
    }
    else if (lo_radius_ < hi_radius_)
    {
        // Cone opens upward
        origin_pos = hi_ - hi_radius_ * dz / dr;
    }

    // Build the cone surface along the given axis
    Real3 origin;
    origin[axis_] = origin_pos;
    build.cone(axis_, neg, origin, dr / dz);

    // Build the bottom and top planes
    if (lo_radius_ >= 0)
    {
        build.plane(axis_, pos, lo_);
    }
    if (hi_radius_ >= 0)
    {
        build.plane(axis_, neg, hi_);
    }

    // Set bounding box
    real_type max_radius = std::max(lo_radius_, hi_radius_);
    for (const Axis& axis : {Axis::x, Axis::y, Axis::z})
    {
        if (axis == axis_)
            continue;
        build.clip_bbox_lower(axis, -max_radius);
        build.clip_bbox_upper(axis, max_radius);
    }
}

//---------------------------------------------------------------------------//
} // namespace celeritas
