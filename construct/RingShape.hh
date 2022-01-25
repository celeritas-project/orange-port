//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RingShape.hh
 * \brief RingShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Cylindrical shell.
 */
class RingShape final : public Shape
{
    using Base = Shape;

  public:
    // Construct ring on the Z axis
    RingShape(real_type inner, real_type outer, real_type lo, real_type hi);

    // Construct ring on the Z axis centered at 0,0,0
    RingShape(real_type inner, real_type outer, real_type height);

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

    // Inner and outer radii
    real_type inner_;
    real_type outer_;

    // Extents along the z axis
    real_type lo_;
    real_type hi_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
