//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/TriangularPrismShape.hh
 * \brief TriangularPrismShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * A triangular prism centered about the origin, cut at the +XY surface
 *
 * The "triangular" side of the prism has a normal pointing into the +x,+y,z=0
 * direction.. This side lies on the z=0 axis. The bounding box of this shape
 * is equivalent to the bounding box of a cuboid with the same x/y/z length
 * inputs.
 */
class TriangularPrismShape final : public Shape
{
    using Base = Shape;

  public:
    // Construct the prism
    TriangularPrismShape(real_type x_length,
                         real_type y_length,
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

    Real3 halfwidths_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
