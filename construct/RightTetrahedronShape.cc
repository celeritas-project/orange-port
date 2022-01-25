//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RightTetrahedronShape.cc
 * \brief RightTetrahedronShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RightTetrahedronShape.hh"

#include "base/Assert.hh"
#include "base/Definitions.hh"
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

using Axis::x;
using Axis::y;
using Axis::z;

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
RightTetrahedronShape::RightTetrahedronShape(real_type x_length,
                                             real_type y_length,
                                             real_type z_length)
    : widths_(x_length, y_length, z_length)
{
    CELER_VALIDATE(x_length > 0.0,
                   << "Triprism x length " << x_length << " must be positive");
    CELER_VALIDATE(y_length > 0.0,
                   << "Triprism y length " << y_length << " must be positive");
    CELER_VALIDATE(z_length > 0.0,
                   << "Triprism z length " << z_length << " must be positive");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string RightTetrahedronShape::type() const
{
    return "right_tetrahedron";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool RightTetrahedronShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief RightTetrahedronShape volume
 */
real_type RightTetrahedronShape::volume() const
{
    real_type vol = widths_[X] * widths_[Y] * widths_[Z] / 6;

    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type RightTetrahedronShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void RightTetrahedronShape::build(ShapeBuilder& build) const
{
    // Build -X, -Y, -Z planes
    build.plane(X, pos, 0.0);
    build.plane(Y, pos, 0.0);
    build.plane(Z, pos, 0.0);

    // Build oblique surface: AB x BC = (-a,b,0) x (0,-b,c) = (bc,ac,ab)
    build.plane(Real3(widths_[Y] * widths_[Z],
                      widths_[X] * widths_[Z],
                      widths_[X] * widths_[Y]),
                Real3(widths_[X], 0, 0));

    // Set upper side of bounding box
    build.clip_bbox_upper(X, widths_[X]);
    build.clip_bbox_upper(Y, widths_[Y]);
    build.clip_bbox_upper(Z, widths_[Z]);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
