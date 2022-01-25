//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/UnionShape.cc
 * \brief UnionShape class definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "UnionShape.hh"

#include <algorithm>
#include <utility>
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"
#include "PlacedShape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
UnionShape::UnionShape(RegionVec interior) : interior_(std::move(interior))
{
    CELER_EXPECT(!interior_.empty());
    CELER_EXPECT(std::all_of(
        interior_.begin(), interior_.end(), [](const Halfspace& hs) {
            return static_cast<bool>(hs.second);
        }));
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string UnionShape::type() const
{
    return "union";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool UnionShape::is_convex() const
{
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Union shape volume
 */
real_type UnionShape::volume() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type UnionShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void UnionShape::build(ShapeBuilder& build_shape) const
{
    // Add an extra push/pop layer in order to create an 'or' node.. The caller
    // will make an 'AND' with a single daughter.
    build_shape.push(CSGCell::LOGIC_OR, {});

    for (const auto& sense_shape : interior_)
    {
        // Apply the local shape's transformation
        build_shape.push(CSGCell::LOGIC_AND, sense_shape.second->transform());

        // Tell unplaced shape to build surfaces
        sense_shape.second->shape()->build(build_shape);

        // Each shape 'AND's the surfaces it builds
        build_shape.pop(sense_shape.first);
    }

    build_shape.pop(inside);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
