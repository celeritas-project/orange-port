//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/IntersectionShape.cc
 * \brief IntersectionShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "IntersectionShape.hh"

#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"
#include "PlacedShape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct from an interior region definition.
 */
IntersectionShape::IntersectionShape(RegionVec interior)
    : interior_(std::move(interior))
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
std::string IntersectionShape::type() const
{
    return "intersection";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool IntersectionShape::is_convex() const
{
    return std::all_of(
        interior_.begin(), interior_.end(), [](const Halfspace& hs) {
            return hs.first == inside && hs.second->shape()->is_convex();
        });
}

//---------------------------------------------------------------------------//
/*!
 * IntersectionShape volume
 */
real_type IntersectionShape::volume() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type IntersectionShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void IntersectionShape::build(ShapeBuilder& build_shape) const
{
    for (const auto& sense_shape : interior_)
    {
        // Apply the local shape's transformation
        build_shape.push(CSGCell::LOGIC_AND, sense_shape.second->transform());

        // Tell unplaced shape to build surfaces
        sense_shape.second->shape()->build(build_shape);

        // Shapes 'AND' together all the surfaces they build
        build_shape.pop(sense_shape.first);
    }
}

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Construct a placed shape (possibly simplified!) from a transformed RDV.
 *
 * An RDV is a "region definition vector": a list of shape + sense
 * intersections defining the interior of an object.
 *
 * The result's embedded shape may be an IntersectionShape *or* the original
 * object.
 */
IntersectionShape::SPConstShape
shape_from_rdv(IntersectionShape::RegionVec interior,
               Transform                    transform,
               ObjectMetadata               md)
{
    CELER_EXPECT(!interior.empty());
    CELER_EXPECT(md);

    PlacedShape::Params params;
    if (interior.size() == 1 && interior.front().first == neg)
    {
        // Copy the underlying shape but apply the transform
        const PlacedShape& daughter = *interior.front().second;
        params.shape                = daughter.shape();
        params.transform            = daughter.transform();
        params.transform.transform(transform);
    }
    else
    {
        // Create intersection shape from RDV
        params.shape = std::make_shared<IntersectionShape>(std::move(interior));
        params.transform = transform;
    }
    params.md = std::move(md);
    return std::make_shared<PlacedShape>(std::move(params));
}

//---------------------------------------------------------------------------//
} // namespace celeritas
