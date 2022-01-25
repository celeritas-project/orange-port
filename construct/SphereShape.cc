//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/SphereShape.cc
 * \brief SphereShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SphereShape.hh"

#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct the SphereShape
 */
SphereShape::SphereShape(real_type radius) : radius_(radius)
{
    CELER_VALIDATE(radius_ > 0, << "Invalid sphere radius " << radius_);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Shape class name
 */
std::string SphereShape::type() const
{
    return "sphere";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool SphereShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sphere volume
 */
real_type SphereShape::volume() const
{
    return sphere_volume(radius_);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Largest sphere radius that fits in this shape
 */
real_type SphereShape::inradius() const
{
    return radius_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add the surfaces
 */
void SphereShape::build(ShapeBuilder& build) const
{
    build.sphere(neg, Real3(0, 0, 0), radius_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
