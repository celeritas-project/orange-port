//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/HopperShape.hh
 * \brief HopperShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * A truncated square pyramid (frustrum)
 *
 * A shape whose top and bottom faces are rectangular parallelepipeds centered
 * about the Z-axis and parallel to the X and Y axes.. It is defined by
 * specifying the half-length of the top face along the X-axis, x_hi,
 * the half-length of the top face along the Y-axis, y_hi, the Z coordinate
 * of the top face, z_hi, the half-length of the bottom face along the X-axis,
 * x_lo, the half-length of the bottom face along the Y-axis, y_lo, and the
 * Z coordinate of the bottom face, z_lo.
 */
class HopperShape final : public Shape
{
    using Base = Shape;

  public:
    // Constructors
    HopperShape(real_type x_hi,
                real_type y_hi,
                real_type z_hi,
                real_type x_lo,
                real_type y_lo,
                real_type z_lo);

    //// DERIVED INTERFACE ////

    // Shape class name
    std::string type() const;

    // Whether the shape is convex (no internal surface crossings)
    bool is_convex() const;

    // Shape volume
    real_type volume() const;

    // Largest sphere radius that fits in this shape
    real_type inradius() const;

    // Build the shape's surfaces
    void build(ShapeBuilder& builder) const;

  private:
    //// DATA ////

    // Lower half-widths
    real_type lo_hx_;
    real_type lo_hy_;

    // Upper half-widths
    real_type hi_hx_;
    real_type hi_hy_;

    // Extents along the Z axis
    real_type lo_z_;
    real_type hi_z_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
