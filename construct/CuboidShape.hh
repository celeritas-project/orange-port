//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CuboidShape.hh
 * \brief CuboidShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

#include <memory>

#include "../Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Rectangular paralellepiped shape.
 */
class CuboidShape final : public Shape
{
    using Base = Shape;

  public:
    //@{
    //! Types
    using PairDbl = std::pair<real_type, real_type>;
    //@}

  public:
    // Construct from alternating lo/hi pairs (keno)
    static CuboidShape from_bounds(PairDbl x, PairDbl y, PairDbl z);

    // Construct from bounds
    CuboidShape(Real3 lo, Real3 hi);

    //// ACCESSORS ////

    //! Lower coordinates
    const Real3& lower() const { return lo_; }

    //! Upper coordinates
    const Real3& upper() const { return hi_; }

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
    Real3 lo_;
    Real3 hi_;
};

//---------------------------------------------------------------------------//
/*!
 * Make a cube shape
 *
 * \todo Possibly replace later with a 'cube' shape to make shape types better
 * reflect user input.
 */
std::shared_ptr<const Shape> make_cube(real_type lo, real_type hi);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
