//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RingShape.cc
 * \brief RingShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RingShape.hh"

#include "base/Assert.hh"
#include "base/Constants.hh"
#include "base/Definitions.hh"
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct the ring
 */
RingShape::RingShape(real_type inner, real_type outer, real_type lo, real_type hi)
    : inner_(inner), outer_(outer), lo_(lo), hi_(hi)
{
    CELER_VALIDATE(inner > 0,
                   << "Inner radius " << inner << " for ring must be positive");
    CELER_VALIDATE(inner < outer,
                   << "Inner radius " << inner
                   << " for ring must be less than outer " << outer);
    CELER_VALIDATE(lo < hi,
                   << "Lower extent " << lo
                   << " must be less than upper extent " << hi << " for ring");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct the ring (SWORD constructor)
 *
 * The result is centered on the origin and along the Z axis
 */
RingShape::RingShape(real_type inner, real_type outer, real_type height)
    : inner_(inner), outer_(outer), lo_(-height / 2), hi_(height / 2)
{
    CELER_VALIDATE(inner > 0,
                   << "Inner radius " << inner
                   << " for cylshell must be positive");
    CELER_VALIDATE(inner < outer,
                   << "Inner radius " << inner
                   << " for cylshell must be less than outer " << outer);
    CELER_VALIDATE(height,
                   << "Height " << height << " must be positive for cylshell");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Shape class name
 */
std::string RingShape::type() const
{
    return "ring";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool RingShape::is_convex() const
{
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief RingShape volume
 */
real_type RingShape::volume() const
{
    using constants::pi;

    // pi * h * (outer_radius^2 - inner_radius^2)
    real_type vol = pi * (hi_ - lo_) * (ipow<2>(outer_) - ipow<2>(inner_));
    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type RingShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void RingShape::build(ShapeBuilder& build) const
{
    // Build inner and outer cylinder
    build.cyl(Axis::z, pos, inner_);
    build.cyl(Axis::z, neg, outer_);

    // Build bounding planes (above lower, below upper)
    build.plane(Axis::z, pos, lo_);
    build.plane(Axis::z, neg, hi_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
