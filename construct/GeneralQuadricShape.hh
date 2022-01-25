//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/GeneralQuadricShape.hh
 * \brief GeneralQuadricShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * A general shape as described by a 10 coefficient quadratic equation
 *
 * is a surface using a quadratic equation of the form:
 * \f[
    aX^2 + bY^2 + cZ^2 + dXY + eYZ + fXZ + gX + hY + iZ + j = 0
    \f]
 *
 * Notice that the e and f coefficients may be switched from typical text book
 * notation.
 */
class GeneralQuadricShape final : public Shape
{
    using Base = Shape;

  public:
    //! Construct the shape
    GeneralQuadricShape(real_type a,
                        real_type b,
                        real_type c,
                        real_type d,
                        real_type e,
                        real_type f,
                        real_type g,
                        real_type h,
                        real_type i,
                        real_type j);

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

    // Second-order terms (a, b, c)
    real_type a_, b_, c_;
    // Second-order cross terms (d, e, f)
    real_type d_, e_, f_;
    // First-order terms (g, h, i)
    real_type g_, h_, i_;
    // Constant term
    real_type j_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
