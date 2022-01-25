//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/PlaneShape.cc
 * \brief PlaneShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "PlaneShape.hh"

#include "base/GeometryUtils.hh"
#include "base/VectorFunctions.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct from normal and coincident point.
 *
 * The normal does not have to be normalized. Its normalization does *not*
 * affect the \c point argument.
 */
PlaneShape::PlaneShape(Real3 normal, Real3 point)
    : normal_(normal), point_(point)
{
    Insist(normal_ != Real3(0, 0, 0), "Invalid null plane direction vector");
}

//---------------------------------------------------------------------------//
/*!
 * Construct the PlaneShape from ax + by + cz - d < 0;
 *
 * KENO specifies "inside" ("positive" sense) as having "az + by + cz + d > 0",
 * whereas "inside" (negative ) in GG/ORANGE is "ax + by + cz - d < 0".. So
 * KENO input must flip the sign of the normal to be consistent with this
 * constructor.
 */
PlaneShape PlaneShape::from_displacement(Real3 normal, real_type displacement)
{
    Real3 point = normal * (displacement / dot_product(normal, normal));
    return PlaneShape(normal, point);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string PlaneShape::type() const
{
    return "plane";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool PlaneShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief PlaneShape volume
 */
real_type PlaneShape::volume() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type PlaneShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void PlaneShape::build(ShapeBuilder& build) const
{
    build.plane(normal_, point_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
