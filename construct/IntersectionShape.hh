//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/IntersectionShape.hh
 * \brief IntersectionShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
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
class ObjectMetadata;
class PlacedShape;
//---------------------------------------------------------------------------//
/*!
 * The "shape" defined by the intersection of other shapes.
 *
 * The intersection shape is used when building holes and boundaries: it is
 * essentially a "region definition vector" as a shape. This allows daughter
 * universe boundaries to be excluded from other regions in the parent.
 *
 * It is also used to define shapes that have "chords" applied to them.
 *
 * \note This is the renamed "MetaShape" from GG (shape of shapes).
 */
class IntersectionShape final : public Shape
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
    explicit IntersectionShape(RegionVec interior);

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
// Construct a placed shape from a transformed RDV
IntersectionShape::SPConstShape
shape_from_rdv(IntersectionShape::RegionVec interior,
               Transform                    transform,
               ObjectMetadata               md);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
