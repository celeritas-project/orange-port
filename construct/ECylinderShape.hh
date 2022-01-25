//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ECylinderShape.hh
 * \brief ECylinderShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

#include "base/Array.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
  * Cylinder along the Z axis with an elliptical cross section.
 *
 * The equation of an Z-aligned elliptical cylinder at the origin is \f[
  \frac{x^2}{r_x^2} + \frac{y^2}{r_y^2} - 1 = 0
  \f]
 */
class ECylinderShape final : public Shape
{
    using Base = Shape;

  public:
    //@{
    //! Typedefs
    using PlaneVector = Array<real_type, 2>;
    //@}

  private:
    //// DATA ////

    // Cylinder x and y radii
    PlaneVector radii_;

    // Extents along the Z axis
    real_type lo_;
    real_type hi_;

  public:
    // Construct the ECylinder
    ECylinderShape(PlaneVector radii, real_type lo_z, real_type hi_z);

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
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
