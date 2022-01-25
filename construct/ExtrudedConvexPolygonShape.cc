//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ExtrudedConvexPolygonShape.cc
 * \brief ExtrudedConvexPolygonShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ExtrudedConvexPolygonShape.hh"

#include <cmath>

#include "base/Assert.hh"
#include "base/Constants.hh"
#include "base/Definitions.hh"
#include "base/Range.hh"
#include "base/GeometryUtils.hh"
#include "detail/PolygonUtils.hh"
#include "detail/ShapeBuilder.hh"

using Axis::x;
using Axis::y;
using Axis::z;

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct the extruded convex polygon
 *
 * Specification of vertices follows the counter clockwise convention
 *
 * Repeating the first coordinate pair as the last coordinate pair helps
 * with reducing complexity of the polygon related methods
 */
ExtrudedConvexPolygonShape::ExtrudedConvexPolygonShape(VecPoint  vertices,
                                                       real_type lo,
                                                       real_type hi)
    : vertices_(std::move(vertices)), lo_(lo), hi_(hi)
{
    CELER_VALIDATE(this->num_polygon_vertices() >= 3,
                   << "Invalid number of polygon vertices "
                   << num_polygon_vertices());
    CELER_VALIDATE(lo <= hi,
                   << "Lower extent " << lo
                   << " must be less than upper extent " << hi
                   << " for extruded convex polygon");
    CELER_VALIDATE(this->is_simple(),
                   << "Vertices must create a simple polygon");
    CELER_VALIDATE(this->vertices_counter_clockwise(),
                   << "Vertices must be specified in a counter-clockwise "
                      "order");
    CELER_VALIDATE(this->is_convex(),
                   << "Vertices must create a convex polygon");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Shape class name
 */
std::string ExtrudedConvexPolygonShape::type() const
{
    return "extruded_convex_polygon";
}

//---------------------------------------------------------------------------//
/*!
 * \brief ExtrudedConvexPolygonShape volume
 */
real_type ExtrudedConvexPolygonShape::volume() const
{
    return polygon_area() * (hi_ - lo_);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Largest sphere radius that fits in this shape
 */
real_type ExtrudedConvexPolygonShape::inradius() const
{
    // find the extents of the polygon
    real_type min_x = vertices_[0][X];
    real_type max_x = vertices_[0][X];
    real_type min_y = vertices_[0][Y];
    real_type max_y = vertices_[0][Y];

    // find max and min x,y extruded convex polygonatic sides
    for (auto n : range(num_polygon_vertices()))
    {
        size_type m = next_vertex(n);

        const Point& vert_n = vertices_[n];
        const Point& vert_m = vertices_[m];

        min_x = std::min(min_x, vert_n[X]);
        max_x = std::max(max_x, vert_n[X]);
        min_y = std::min(min_y, vert_n[Y]);
        max_y = std::max(max_y, vert_n[Y]);
    }

    real_type hmax    = half * std::max(max_x - min_x, max_y - min_y);
    real_type hheight = half * (hi_ - lo_);
    return std::sqrt(ipow<2>(hheight) + 2.0 * hmax);
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void ExtrudedConvexPolygonShape::build(ShapeBuilder& build) const
{
    real_type min_x = vertices_[0][X];
    real_type max_x = vertices_[0][X];
    real_type min_y = vertices_[0][Y];
    real_type max_y = vertices_[0][Y];

    // Build extruded convex polygonatic sides
    for (auto n : range(num_polygon_vertices()))
    {
        size_type m = next_vertex(n);

        const Point& vert_n = vertices_[n];
        const Point& vert_m = vertices_[m];

        Real3 point(vert_n[X], vert_n[Y], 0.0);
        Real3 normal(-(vert_n[Y] - vert_m[Y]), (vert_n[X] - vert_m[X]), 0.0);
        build.plane(normal, point);

        min_x = std::min(min_x, vert_n[X]);
        max_x = std::max(max_x, vert_n[X]);
        min_y = std::min(min_y, vert_n[Y]);
        max_y = std::max(max_y, vert_n[Y]);
    }

    // Build top and bottom
    build.plane(Z, pos, lo_);
    build.plane(Z, neg, hi_);

    // Set the +/- x & y axis bounding box
    build.clip_bbox_lower(X, min_x);
    build.clip_bbox_upper(X, max_x);
    build.clip_bbox_lower(Y, min_y);
    build.clip_bbox_upper(Y, max_y);
}

//---------------------------------------------------------------------------//
/*!
 * Index of next vertex in sequence after vertex n
 */
size_type ExtrudedConvexPolygonShape::next_vertex(const size_type n) const
{
    CELER_EXPECT(num_polygon_vertices() > 1);
    CELER_EXPECT(n < num_polygon_vertices());
    return (n + 1) % num_polygon_vertices();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Index of previous vertex in sequence after vertex n
 */
size_type ExtrudedConvexPolygonShape::prev_vertex(const size_type n) const
{
    CELER_EXPECT(num_polygon_vertices() > 1);
    CELER_EXPECT(n < num_polygon_vertices());
    // Adding by num_polygon_vertices() because we're dealing with
    // unsigned ints and if n = 0, then n - 1 would cause errors
    return ((n + num_polygon_vertices()) - 1) % num_polygon_vertices();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether any polygon edges overlap
 */
bool ExtrudedConvexPolygonShape::edges_overlap() const
{
    using detail::edges_overlap;

    CELER_EXPECT(num_polygon_vertices() >= 3);

    for (auto n : range(num_polygon_vertices()))
    {
        if (edges_overlap(vertices_[prev_vertex(n)],
                          vertices_[n],
                          vertices_[next_vertex(n)]))
        {
            return true;
        }
    }
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether any polygon edges intersect at points other than vertices
 *
 *  NOTE: this is an O(N^2) operation so it can be expensive
 */
bool ExtrudedConvexPolygonShape::edges_intersect() const
{
    using detail::edges_intersect;

    CELER_EXPECT(num_polygon_vertices() >= 3);

    for (auto n : range(num_polygon_vertices()))
    {
        for (auto m : range(num_polygon_vertices()))
        {
            // Ignore cases considering the same edge
            if (n == m)
                continue;
            // Ignore cases where edges share a vertex
            if (next_vertex(n) == m)
                continue;
            // Ignore cases where edges share a vertex
            if (prev_vertex(n) == m)
                continue;

            if (edges_intersect(vertices_[n],
                                vertices_[next_vertex(n)],
                                vertices_[m],
                                vertices_[next_vertex(m)]))
            {
                return true;
            }
        }
    }
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether polygon is simple
 */
bool ExtrudedConvexPolygonShape::is_simple() const
{
    return !(edges_overlap() || edges_intersect());
}

//---------------------------------------------------------------------------//
/*!
 * \brief The (signed) area of the defining polygon
 *
 * For details: http://mathworld.wolfram.com/PolygonArea.html
 */
real_type ExtrudedConvexPolygonShape::polygon_area() const
{
    real_type area = 0;
    for (auto n : range(num_polygon_vertices()))
    {
        area += vertices_[n][Axis::x] * vertices_[next_vertex(n)][Axis::y]
                - vertices_[n][Axis::y] * vertices_[next_vertex(n)][Axis::x];
    }
    CELER_ASSERT(std::fabs(area) > 1.0e-8);
    return half * area;
}

//---------------------------------------------------------------------------//
/*!
 * Whether polygon vertices are arranged in a Counter Clockwise order
 */
bool ExtrudedConvexPolygonShape::vertices_counter_clockwise() const
{
    CELER_EXPECT(is_simple());
    return polygon_area() > 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Whether polygon is convex
 *
 * Evaluation loosely based on Graham scan method wherein any three points
 * constituting a right turn implies a non-convex polygon.
 *
 * For details: https://en.wikipedia.org/wiki/Graham_scan
 */
bool ExtrudedConvexPolygonShape::is_convex() const
{
    using detail::is_right_turn;

    CELER_EXPECT(vertices_counter_clockwise());
    for (auto n : range(num_polygon_vertices()))
    {
        if (is_right_turn(vertices_[prev_vertex(n)],
                          vertices_[n],
                          vertices_[next_vertex(n)]))
            return false;
    }
    return true;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
