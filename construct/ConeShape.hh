//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ConeShape.hh
 * \brief ConeShape class declaration
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
 * Truncated cone shape.
 *
 * The cone shape is aligned along the given axis; it has lower and upper
 * radii. Either of these can be set to zero to make a non-truncated
 * single-sheet cone.
 */
class ConeShape final : public Shape
{
    using Base = Shape;

  public:
    //@{
    //! Type aliases
    using PairDbl = std::pair<real_type, real_type>;
    //@}

  public:
    // Construct the cone shape
    ConeShape(def::XYZ axis, PairDbl radii, PairDbl extents);

    //// ACCESSORS ////

    //! High radius
    real_type high_radius() const { return hi_radius_; }

    //! Low radius
    real_type low_radius() const { return lo_radius_; }

    //! High extent along z-axis
    real_type high_extent() const { return hi_; }

    //! Low extent along z-axis
    real_type low_extent() const { return lo_; }

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
    // Axis of the cone
    def::XYZ axis_;

    // Radii
    real_type lo_radius_;
    real_type hi_radius_;

    // Extents along the given axis
    real_type lo_;
    real_type hi_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
