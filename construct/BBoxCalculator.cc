//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/BBoxCalculator.cc
 * \brief BBoxCalculator class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "BBoxCalculator.hh"

#include "orange/surfaces/SurfaceContainer.hh"
#include "PlacedShape.hh"
#include "Shape.hh"
#include "CSGTree.hh"
#include "detail/ShapeBuilder.hh"
#include "detail/SurfaceInserter.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Calculate the bounding box for a transformed shape.
 *
 * This reuses the same surface container as previous calls to the operator.
 */
auto BBoxCalculator::operator()(const Shape& s, const Transform& t)
    -> result_type
{
    SurfaceContainer        surfaces;
    detail::SurfaceInserter insert_surface(&surfaces);
    CSGTree                 tree;
    detail::ShapeBuilder    build_shape(insert_surface, tree);

    // Add the user-provided shape: inside all given surfaces
    build_shape.push(CSGCell::LOGIC_AND, t);
    s.build(build_shape);
    build_shape.pop(neg);

    // Save result
    auto result = build_shape();
    return result.bbox;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the bounding box for a shape.
 */
auto BBoxCalculator::operator()(const Shape& s) -> result_type
{
    return (*this)(s, Transform{});
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the bounding box for a placed shape.
 */
auto BBoxCalculator::operator()(const PlacedShape& ps) -> result_type
{
    return (*this)(*ps.shape(), ps.transform());
}

//---------------------------------------------------------------------------//
} // namespace celeritas
