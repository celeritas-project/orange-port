//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/WedgeShape.cc
 * \brief WedgeShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "WedgeShape.hh"

#include <cmath>
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
WedgeShape::WedgeShape(real_type x_length,
                       real_type x_base_pt,
                       real_type y_base_pt,
                       real_type z_length)
    : x_length_(x_length)
    , z_length_(z_length)
    , x_pt_(x_base_pt)
    , y_pt_(y_base_pt)
{
    CELER_VALIDATE(x_length_ > 0.0,
                   << "Non-positive x-length in wedge: " << x_length_);
    CELER_VALIDATE(z_length_ > 0.0,
                   << "Non-positive z-length in wedge: " << z_length_);
    CELER_VALIDATE(x_pt_ >= 0.0 && x_pt_ <= x_length_,
                   << "Base point in wedge be somewhere along the x-length: "
                   << x_pt_);
    CELER_VALIDATE(y_pt_ > 0.0,
                   << "Base point in wedge must be greater than zero: "
                   << y_pt_);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string WedgeShape::type() const
{
    return "wedge";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool WedgeShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief WedgeShape volume
 */
real_type WedgeShape::volume() const
{
    real_type vol = half * x_length_ * z_length_ * y_pt_;
    CELER_ENSURE(vol > 0.0);

    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type WedgeShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void WedgeShape::build(ShapeBuilder& build) const
{
    // Build side surface along x-axis
    build.plane(Axis::y, pos, 0.0);

    // Build first oblique surface
    build.plane(Real3(-y_pt_, x_pt_, 0.0), Real3(0.0, 0.0, 0.0));

    // Build second oblique surface
    build.plane(Real3(y_pt_, x_length_ - x_pt_, 0.0),
                Real3(x_length_, 0.0, 0.0));

    // Build bottom surface.
    build.plane(Axis::z, pos, 0.0);

    // Build top surface
    build.plane(Axis::z, neg, z_length_);

    // Set the +/- x axis bounding box
    build.clip_bbox_lower(Axis::x, 0.0);
    build.clip_bbox_upper(Axis::x, x_length_);
    // Set the +y axis bounding box
    build.clip_bbox_upper(Axis::y, y_pt_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
