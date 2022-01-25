//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/FaceNameCalculator.cc
 * \brief FaceNameCalculator class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "FaceNameCalculator.hh"

#include "orange/surfaces/SurfaceContainer.hh"
#include "Shape.hh"
#include "CSGTree.hh"
#include "detail/SurfaceInserter.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Determine the names of the faces of the given shape.
 */
auto FaceNameCalculator::operator()(const Shape& shape) -> result_type
{
    // Construct shape builder
    SurfaceContainer        surfaces;
    detail::SurfaceInserter insert_surface(&surfaces);
    CSGTree                 tree;
    detail::ShapeBuilder    build_shape(insert_surface, tree);

    // Add callback for saving the surface name
    std::vector<std::string> surface_names;
    build_shape.set_surface_callback(
        [&surface_names](SurfaceId, std::string name) {
            surface_names.push_back(std::move(name));
        });

    // Add the shape to construct the surface names
    build_shape.push(CSGCell::LOGIC_AND, Transform{});
    shape.build(build_shape);
    build_shape.pop();

    CELER_ENSURE(!surface_names.empty());
    return surface_names;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
