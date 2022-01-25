//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/UnionShape.hh
 * \brief UnionShape class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

#include <memory>
#include <utility>
#include <vector>

#include "orange/Definitions.hh"

namespace celeritas
{
class PlacedShape;
//---------------------------------------------------------------------------//
/*!
 * The "shape" defined by the union of other shapes.
 */
class UnionShape final : public Shape
{
    using Base = Shape;

  public:
    //@{
    //! Public type aliases
    using SPConstShape = std::shared_ptr<const PlacedShape>;
    using Halfspace    = std::pair<Sense, SPConstShape>;
    using RegionVec    = std::vector<Halfspace>;
    //@}

  public:
    // Constructors
    explicit UnionShape(RegionVec interior);

    //// ACCESSORS ////

    //! Access the region definition vector defining this shape.
    const RegionVec& interior() const { return interior_; }

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
    RegionVec interior_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
