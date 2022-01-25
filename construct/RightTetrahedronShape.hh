//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RightTetrahedronShape.hh
 * \brief RightTetrahedronShape class declaration
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
 * A tetrahedron in the +XYZ octant with three edges along X, Y, and Z
 */
class RightTetrahedronShape final : public Shape
{
    using Base = Shape;

  public:
    // Construct the tet
    RightTetrahedronShape(real_type x_length,
                          real_type y_length,
                          real_type z_length);

    //// ACCESSORS ////

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

    Real3 widths_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
