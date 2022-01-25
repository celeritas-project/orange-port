//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ChordShape.hh
 * \brief ChordShape class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Single axis-aligned plane.
 */
class ChordShape final : public Shape
{
    using Base = Shape;

  public:
    // Constructor
    ChordShape(def::XYZ axis, Sense sense, real_type position);

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
    void build(ShapeBuilder& build) const;

  private:
    //// DATA ////

    def::XYZ  axis_;
    Sense     sense_;
    real_type position_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
