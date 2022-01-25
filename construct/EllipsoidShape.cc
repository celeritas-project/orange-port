//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/EllipsoidShape.cc
 * \brief EllipsoidShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "EllipsoidShape.hh"

#include <algorithm>

#include "base/Assert.hh"
#include "base/Constants.hh"
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
/*!
 * Constructor.
 */
EllipsoidShape::EllipsoidShape(Real3 radii) : radii_(radii)
{
    for (Axis ax : {X, Y, Z})
    {
        CELER_VALIDATE(radii_[ax] >= 0,
                       << "Ellipsoid axis along XYZ"[ax]
                       << " must be positive (given r=" << radii_[ax] << ")");
    }
}
//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string EllipsoidShape::type() const
{
    return "ellipsoid";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool EllipsoidShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Shape volume
 */
real_type EllipsoidShape::volume() const
{
    using constants::four_pi;
    real_type vol = four_pi / 3.0 * radii_[X] * radii_[Y] * radii_[Z];
    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type EllipsoidShape::inradius() const
{
    real_type max_radius = std::max(radii_[X], radii_[Y]);
    return std::max(max_radius, radii_[Z]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct surfaces for this shape
 */
void EllipsoidShape::build(ShapeBuilder& build) const
{
    // Second-order coefficients are product of the other two squared radii;
    // Zeroth-order coefficient is the product of all three squared radii
    Real3 rsq;
    for (Axis ax : {X, Y, Z})
    {
        rsq[ax] = radii_[ax] * radii_[ax];
    }

    Real3     abc(1, 1, 1);
    real_type g = -1;
    for (Axis ax : {X, Y, Z})
    {
        g *= rsq[ax];
        for (Axis nax : {X, Y, Z})
        {
            if (ax != nax)
            {
                abc[ax] *= rsq[nax];
            }
        }
    }

    // Build the quadric
    build.simple_quadric(abc, Real3(0, 0, 0), g, Real3(0, 0, 0));

    // Clip radial directions based on outer shape radius
    for (Axis ax : {X, Y, Z})
    {
        build.clip_bbox_lower(ax, -radii_[ax]);
        build.clip_bbox_upper(ax, radii_[ax]);
    }
}

//---------------------------------------------------------------------------//
} // namespace celeritas
