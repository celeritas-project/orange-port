//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ChordShape.cc
 * \brief ChordShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ChordShape.hh"

#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct the ChordShape
 */
ChordShape::ChordShape(def::XYZ axis, Sense sense, real_type position)
    : axis_(axis), sense_(sense), position_(position)
{
    CELER_VALIDATE(axis < def::END_XYZ, << "Invalid axis " << axis);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Shape name
 */
std::string ChordShape::type() const
{
    return "chord";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool ChordShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief PlaneShape volume
 */
real_type ChordShape::volume() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type ChordShape::inradius() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
void ChordShape::build(ShapeBuilder& build) const
{
    // Build surfaces
    build.plane(axis_, sense_, position_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
