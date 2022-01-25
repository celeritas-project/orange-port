//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/Shape.hh
 * \brief Shape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <string>

namespace celeritas
{
namespace detail
{
class ShapeBuilder;
}

//---------------------------------------------------------------------------//
/*!
 * Base class for a 'factory' that defines and can build a shape.
 *
 * Shapes are the region of space defined by a primitive object in a KENO
 * input, a collection of intersected surfaces. These shape classes are *only*
 * used for construction and testing. The internal representation of cells will
 * only consist of surfaces.
 *
 * The "inradius" (http://mathworld.wolfram.com/Insphere.html) is used for
 * particle/pebble sampling: it's the largest sphere that can fit inside this
 * shape.
 *
 * \todo Use sentinel instead of zero for cell/inradius.
 */
class Shape
{
  public:
    //@{
    //! Types
    using ShapeBuilder = detail::ShapeBuilder;
    //@}

  public:
    //! Virtual destructor
    virtual ~Shape() = 0;

    //// ACCESSORS ////

    //! Description of the type of shape
    virtual std::string type() const = 0;

    //! Whether the shape is convex (no internal surface crossings)
    virtual bool is_convex() const = 0;

    //! Shape's interior cell [cm^3]
    virtual real_type volume() const = 0;

    //! Radius of the largest sphere that can fit in this shape [cm]
    virtual real_type inradius() const = 0;

    //// CONSTRUCTION ////

    // Construct surfaces in the given unit
    virtual void build(ShapeBuilder& builder) const = 0;

    // TODO: add intersection tests for improved construction
};

//---------------------------------------------------------------------------//
// Print a shape
std::ostream& operator<<(std::ostream&, const Shape&);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
