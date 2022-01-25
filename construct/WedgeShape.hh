//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/WedgeShape.hh
 * \brief WedgeShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Wedge shape is a right-triangular prism having five faces
 *
 * The two ends are triangles, and the three sides are rectangles.
 *
 * One side is in the XZ plane at Y = 0, and the bottom face is in the
 * XY plane at Z = 0, with a corner at the origin.
 *
 * It is defined by specifying the length of the base along the X-axis,
 * x_length, the X and Y coordinate where the other two sides meet, x_base_pt
 * and y_base_pt, and the length along the Z-axis, z_length.
 *
 * It is restricted to the first quadrant of the XY plane.
 *
 * I.e., x_base_pt must be 0 or greater and y_base_pt must be strictly greater
 * than 0.
 */
class WedgeShape final : public Shape
{
    using Base = Shape;

  public:
    // Constructors
    WedgeShape(real_type x_length,
               real_type x_base_pt,
               real_type y_base_pt,
               real_type z_length);

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

    // Length along x- and z-axis
    real_type x_length_;
    real_type z_length_;

    // (X,Y) corner point in wedge base
    real_type x_pt_;
    real_type y_pt_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
