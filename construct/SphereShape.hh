//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/SphereShape.hh
 * \brief SphereShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * A simple spherical volume.
 */
class SphereShape final : public Shape
{
    using Base = Shape;

  public:
    // Construct the sphere
    explicit SphereShape(real_type radius);

    //// ACCESSORS ////

    //! Sphere radius
    real_type radius() const { return radius_; }

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

    real_type radius_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
