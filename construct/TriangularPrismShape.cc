//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/TriangularPrismShape.cc
 * \brief TriangularPrismShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "TriangularPrismShape.hh"

#include "orange/Definitions.hh"
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
TriangularPrismShape::TriangularPrismShape(real_type x_length,
                                           real_type y_length,
                                           real_type z_length)
    : halfwidths_(x_length / 2, y_length / 2, z_length / 2)
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
std::string TriangularPrismShape::type() const
{
    return "triangular_prism";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool TriangularPrismShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief TriangularPrismShape volume
 */
real_type TriangularPrismShape::volume() const
{
    real_type vol = 4 * halfwidths_[X] * halfwidths_[Y] * halfwidths_[Z];

    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type TriangularPrismShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void TriangularPrismShape::build(ShapeBuilder& build) const
{
    // Build -X and -Y planes
    build.plane(X, pos, -halfwidths_[X]);
    build.plane(Y, pos, -halfwidths_[Y]);

    // Build oblique surface at the origin
    Real3 normal(halfwidths_[Y], halfwidths_[X], 0);
    build.plane(normal, Real3(0, 0, 0));

    // Build top and bottom surfaces
    build.plane(Z, pos, -halfwidths_[Z]);
    build.plane(Z, neg, halfwidths_[Z]);

    // Set upper side of bounding box
    build.clip_bbox_upper(X, halfwidths_[X]);
    build.clip_bbox_upper(Y, halfwidths_[Y]);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
