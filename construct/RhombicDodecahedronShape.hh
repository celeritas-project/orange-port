//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RhombicDodecahedronShape.hh
 * \brief RhombicDodecahedronShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * A regular polyhedron whose surface consists of 12 equal rhombii.
 *
 * The shape is defined by an inscribing radius. This shape is the principal
 * component of a dodecahedral array.
 */
class RhombicDodecahedronShape final : public Shape
{
    using Base = Shape;

  public:
    // Construct the Dodecahedron
    explicit RhombicDodecahedronShape(real_type radius);

    //// ACCESSORS ////

    //! Apothem (inner radius)
    real_type apothem() const { return radius_; }

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

    // the inscribing radius (apothem) of the rhombic dodecahedron
    real_type radius_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
