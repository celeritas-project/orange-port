//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/SlabShape.hh
 * \brief SlabShape class declaration
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
 * Inifinite slab with extents along some axis.
 */
class SlabShape final : public Shape
{
    using Base = Shape;

  public:
    // Construct the SlabShape
    SlabShape(def::XYZ axis, real_type lo, real_type hi);

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

    // Axis normal to the slab faces
    def::XYZ axis_;

    // Extents along the given axis
    real_type lo_;
    real_type hi_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
