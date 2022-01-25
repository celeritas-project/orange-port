//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CylinderShape.hh
 * \brief CylinderShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

#include "base/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Axis-aligned cylinder with finite extents.
 */
class CylinderShape final : public Shape
{
    using Base = Shape;

  public:
    // Constructors
    // Construct the Cylinder
    CylinderShape(def::XYZ axis, real_type radius, real_type lo, real_type hi);

    // Construct a z-aligned cylinder centered on the origin
    CylinderShape(real_type radius, real_type height);

    //// ACCESSORS ////

    //! Axis along the cylinder center
    def::XYZ axis() const { return axis_; }

    //! Cylinder radius
    real_type radius() const { return radius_; }

    //! Coordinate of the cylinder's base along the axis
    real_type lo() const { return lo_; }

    //! Coordinate of the cylinder's top along the axis
    real_type hi() const { return hi_; }

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

    // Axis of the cylinder
    def::XYZ axis_;

    // Cylinder radius
    real_type radius_;

    // Extents along the given axis
    real_type lo_;
    real_type hi_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
