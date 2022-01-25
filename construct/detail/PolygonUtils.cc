//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/PolygonUtils.cc
 * \brief Polygon Utility function declarations
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "PolygonUtils.hh"

#include <cmath>
#include "base/SoftEqual.hh"
#include "base/SoftEquivalence.hh"
#include "base/Assert.hh"
#include "base/Range.hh"
#include "orange/Fuzziness.hh"

using Point = Array<real_type, 2>;
using Axis::x;
using Axis::y;
using std::hypot;

namespace celeritas
{
namespace detail
{
namespace
{
//---------------------------------------------------------------------------//
// INTERNAL HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Whether values are nearly equivalent
 */
bool soft_equiv(real_type a, real_type b)
{
    SoftEqual<real_type> is_soft_equiv{fuzziness().shape_enclosure_rel()};
    return is_soft_equiv(a, b);
}

//---------------------------------------------------------------------------//
/*!
 * Whether points are nearly equivalent
 */
bool soft_equiv(const Point& point_A, const Point& point_B)
{
    real_type distance
        = std::hypot(point_A[X] - point_B[X], point_A[Y] - point_B[Y]);
    SoftZero<real_type> is_soft_zero{fuzziness().shape_enclosure_rel()};
    return is_soft_zero(distance);
}
//---------------------------------------------------------------------------//
} // namespace

//---------------------------------------------------------------------------//
/*!
 * Whether edges associated with the points overlap
 *
 * For each vertex, compute cross and dot of adjacent edges to make
 * sure the edges are not overlapping
 */
bool edges_overlap(const Point& point_A,
                   const Point& point_B,
                   const Point& point_C)
{
    CELER_EXPECT(!soft_equiv(point_A, point_B));
    CELER_EXPECT(!soft_equiv(point_B, point_C));
    const Point     dU       = point_A - point_B;
    const Point     dV       = point_C - point_B;
    const real_type kDotProd = dU[X] * dV[X] + dU[Y] * dV[Y];
    return (kDotProd > 0 && soft_equiv(dU[X] * dV[Y], dU[Y] * dV[X]));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether edges associated with the points intersect
 *
 * For details: https://en.wikipedia.org/wiki/Line-line_intersection
 */
bool edges_intersect(const Point& point_A,
                     const Point& point_B,
                     const Point& point_C,
                     const Point& point_D)
{
    CELER_EXPECT(!soft_equiv(point_A, point_B));
    CELER_EXPECT(!soft_equiv(point_C, point_D));
    // Same line test
    if (soft_equiv(point_A, point_C) && soft_equiv(point_B, point_D))
        return true;
    if (soft_equiv(point_A, point_D) && soft_equiv(point_B, point_C))
        return true;
    const Point     dAB    = point_A - point_B;
    const Point     dAC    = point_A - point_C;
    const Point     dCD    = point_C - point_D;
    const real_type denom1 = dAB[X] * dCD[Y];
    const real_type denom2 = dAB[Y] * dCD[X];
    if (soft_equiv(denom1, denom2))
        return false;
    const real_type denom = denom1 - denom2;
    const real_type t = ((dAC[X]) * (dCD[Y]) - (dAC[Y]) * (dCD[X])) / denom;
    const real_type u = -((dAB[X]) * (dAC[Y]) - (dAB[Y]) * (dAC[X])) / denom;
    // Infinite line Intersection outside a line segment
    if (t < 0.0 || t > 1.0)
        return false;
    if (u < 0.0 || u > 1.0)
        return false;
    // Intersection inside both line segments
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * Whether the three points constitute a right turn
 *
 * This assumes the points are being provided in a counter clockwise manner
 *
 * For details: https://en.wikipedia.org/wiki/Graham_scan
 */
bool is_right_turn(const Point& point_A,
                   const Point& point_B,
                   const Point& point_C)
{
    Point dBA = point_B - point_A;
    Point dCB = point_C - point_B;
    return dBA[X] * dCB[Y] < dBA[Y] * dCB[X]
           && !soft_equiv(dBA[X] * dCB[Y], dBA[Y] * dCB[X]);
}

//---------------------------------------------------------------------------//
} // end namespace detail
} // end namespace celeritas
