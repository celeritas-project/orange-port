//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/HopperShape.cc
 * \brief HopperShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "HopperShape.hh"

#include <algorithm>
#include <cmath>

#include "base/Assert.hh"
#include "base/Constants.hh"
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
HopperShape::HopperShape(real_type x_hi,
                         real_type y_hi,
                         real_type z_hi,
                         real_type x_lo,
                         real_type y_lo,
                         real_type z_lo)
    : lo_hx_(x_lo)
    , lo_hy_(y_lo)
    , hi_hx_(x_hi)
    , hi_hy_(y_hi)
    , lo_z_(z_lo)
    , hi_z_(z_hi)
{
    CELER_VALIDATE(x_lo > 0.0,
                   << "Hopper x_lo " << x_lo << " must be positive");
    CELER_VALIDATE(y_lo > 0.0,
                   << "Hopper y_lo " << y_lo << " must be positive");
    CELER_VALIDATE(x_hi > 0.0,
                   << "Hopper x_hi " << x_hi << " must be positive");
    CELER_VALIDATE(y_hi > 0.0,
                   << "Hopper y_hi " << y_hi << " must be positive");
    CELER_VALIDATE(z_hi > z_lo,
                   << "Hopper z_hi " << z_hi << " must be greater than z_lo"
                   << z_lo);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string HopperShape::type() const
{
    return "hopper";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool HopperShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief HopperShape volume
 */
real_type HopperShape::volume() const
{
    const real_type bottom = lo_hx_ * lo_hy_ * 4;
    const real_type top    = hi_hx_ * hi_hy_ * 4;

    real_type vol = (hi_z_ - lo_z_) / 3
                    * (bottom + std::sqrt(bottom * top) + top);

    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type HopperShape::inradius() const
{
    // Determine max half x and y
    real_type mhx = std::max(lo_hx_, hi_hx_);
    real_type mhy = std::max(lo_hy_, hi_hy_);

    real_type hz = half * (hi_z_ - lo_z_);
    return std::sqrt(ipow<2>(hz) + ipow<2>(mhx) + ipow<2>(mhy));
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void HopperShape::build(ShapeBuilder& build) const
{
    Real3 delta(hi_hx_, hi_hy_, hi_z_);
    delta -= Real3(lo_hx_, lo_hy_, lo_z_);

    Real3 origin(0, 0, lo_z_);

    // +X side
    origin[X] = lo_hx_;
    build.plane(Real3(delta[Z], 0, -delta[X]), origin);

    // +Y side
    origin[Y] = lo_hy_;
    build.plane(Real3(0, delta[Z], -delta[Y]), origin);

    // -X side
    origin[X] = -lo_hx_;
    build.plane(Real3(-delta[Z], 0, -delta[X]), origin);

    // -Y side
    origin[Y] = -lo_hy_;
    build.plane(Real3(0, -delta[Z], -delta[Y]), origin);

    // Build top and bottom surfaces
    build.plane(Z, pos, lo_z_);
    build.plane(Z, neg, hi_z_);

    // Clip xy directions based on hi/lo
    const real_type max_x = std::max(lo_hx_, hi_hx_);
    const real_type max_y = std::max(lo_hy_, hi_hy_);
    build.clip_bbox_lower(X, -max_x);
    build.clip_bbox_upper(X, max_x);
    build.clip_bbox_lower(Y, -max_y);
    build.clip_bbox_upper(Y, max_y);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
