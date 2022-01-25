//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/EllipsoidShape.hh
 * \brief EllipsoidShape class declaration
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
  * Axis-aligned ellipsoid centered at the origin.
 *
 * The equation of an axis-aligned ellipsoid at the origin is \f[
  \frac{x^2}{r_x^2} + \frac{y^2}{r_y^2} + \frac{z^2}{r_z^2} - 1 = 0
  \f]
 */
/*!
 * \example celeritas/construct/test/tstEllipsoidShape.cc
 */
class EllipsoidShape final : public Shape
{
    using Base = Shape;

  private:
    // Radii
    Real3 radii_;

  public:
    // Construct from radii along each axis
    explicit EllipsoidShape(Real3 radii);

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
