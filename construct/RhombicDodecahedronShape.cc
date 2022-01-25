//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RhombicDodecahedronShape.cc
 * \brief RhombicDodecahedronShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RhombicDodecahedronShape.hh"

#include "base/Assert.hh"
#include "base/Constants.hh"
#include "base/Face.hh"
#include "base/Range.hh"
#include "orange/surfaces/Plane.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct from the inner radius of the dodecahedron.
 */
RhombicDodecahedronShape::RhombicDodecahedronShape(real_type radius)
    : radius_(radius)
{
    CELER_VALIDATE(radius_ > 0.0,
                   << "rhombic dodecahedron inscribing radius " << radius_
                   << " must be positive");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string RhombicDodecahedronShape::type() const
{
    return "rhombic_dodecahedron";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool RhombicDodecahedronShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief RhombicDodecahedronShape volume
 */
real_type RhombicDodecahedronShape::volume() const
{
    // vol = 16 * sqrt(3) / 9 * a^3
    // where edge length a = (r * 3) / sqrt(6)
    // where inscribing radius r = a / 3 * sqrt(6)
    // vol = 16 * sqrt(3) / 9 * (r * 3) / sqrt(6)^3
    // simplifies to vol = 4 * sqrt(2) * r^3
    constexpr real_type coeff = 4.0 * constants::sqrt_two;
    real_type           vol   = coeff * ipow<3>(radius_);
    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type RhombicDodecahedronShape::inradius() const
{
    return radius_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct surfaces for this shape.
 */
void RhombicDodecahedronShape::build(ShapeBuilder& build) const
{
    using Face_t = RhombicDodecahedronFace<int>;
    for (auto i : range(Face_t::num_faces()))
    {
        // Build "inside" this face with its given name (replacing +/- with
        // p/m), at the given normal, at a distance of the apothem from the
        // origin.
        Face_t      face{i};
        std::string face_name = to_string(face);
        CELER_ASSERT(face_name.size() == 2);
        face_name.front() = face.is_positive() ? 'p' : 'm';
        build.surface(std::move(face_name),
                      Sense::inside,
                      Plane(face.calc_cartesian_normal(), radius_));
    }

    // Calculate the z bounding box extent
    const real_type z_extent = constants::sqrt_two * radius_;

    // Set upper sides of bounding box
    build.clip_bbox_upper(Axis::x, radius_);
    build.clip_bbox_upper(Axis::y, radius_);
    build.clip_bbox_upper(Axis::z, z_extent);
    // Set lower sides of bounding box
    build.clip_bbox_lower(Axis::x, -radius_);
    build.clip_bbox_lower(Axis::y, -radius_);
    build.clip_bbox_lower(Axis::z, -z_extent);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
