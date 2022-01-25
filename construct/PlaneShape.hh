//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/PlaneShape.hh
 * \brief PlaneShape class declaration
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
 * Single arbitrary plane.
 *
 * "Inside" a plane, for a normalized outward \f$n\f$ and point \f$p\f$,
 * \f[
 * \vec n \vd \vec x - \vec n \vd \vec p < 0
 * \f]
 *
 * this differs in sign of the normal component from the KENO definition, where
 * "inside" is
 * \f[
 * \vec n \vd \vec x + d > 0
 * \f]
 */
class PlaneShape final : public Shape
{
    using Base = Shape;

  public:
    static PlaneShape from_displacement(Real3 normal, real_type displacement);

    // Constructors
    PlaneShape(Real3 normal, Real3 point);

    //// ACCESSORS ////

    //! Outward normal
    const Real3& normal() const { return normal_; }

    //! Intersection point
    const Real3& point() const { return point_; }

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
    Real3 normal_;
    Real3 point_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
