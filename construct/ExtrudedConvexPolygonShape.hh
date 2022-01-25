//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ExtrudedConvexPolygonShape.hh
 * \brief ExtrudedConvexPolygonShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

#include <vector>

#include "base/Definitions.hh"
#include "base/Array.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Generic Convex Polygon extruded into Z.
 */
class ExtrudedConvexPolygonShape final : public Shape
{
    using Base = Shape;
    //@{
    //! Typdefs
    using Point        = Array<real_type, 2>;
    using VecPoint     = std::vector<Point>;
    using PointVecDiff = Array<real_type, 2>;
    //@}
  private:
    //// DATA ////

    // Vertices of convex polygon
    VecPoint vertices_;

    // Extents along the z axis
    real_type lo_;
    real_type hi_;

  public:
    // Constructor
    ExtrudedConvexPolygonShape(VecPoint vertices, real_type lo, real_type hi);

    //// ACCESSORS ////

    //! Vertices of convex polygon
    const VecPoint& vertices() const { return vertices_; }

    //! Number of vertices comprising the extruded polygon
    size_type num_polygon_vertices() const { return vertices_.size(); }

    //! Lower point
    real_type lower() const { return lo_; }

    //! Upper point
    real_type upper() const { return hi_; }

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
    // Index of next vertex in sequence after vertex n
    size_type next_vertex(const size_type) const;

    // Index of previous vertex in sequence after vertex n
    size_type prev_vertex(const size_type) const;

    // Whether any polygon edges overlap
    bool edges_overlap() const;

    // Whether any polygon edges intersect
    bool edges_intersect() const;

    // Whether polygon is simple
    bool is_simple() const;

    // The (signed) area of the defining polygon
    real_type polygon_area() const;

    // Whether polygon vertices are arranged in a Counter Clockwise order
    bool vertices_counter_clockwise() const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
