//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/PrismShape.cc
 * \brief PrismShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "PrismShape.hh"

#include <cmath>
#include <limits>

#include "base/Constants.hh"
#include "base/Range.hh"
#include "base/Definitions.hh"
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

using constants::pi;
using constants::two_pi;

using Axis::x;
using Axis::y;
using Axis::z;

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct the prism
 *
 * The apothem is the distance from center to midpoint of any of its sides.
 * This value is less than the circumradius: \f[

    a = r \cos(\pi / N)

 \f]
 * where \em N is the number of sides.
 */
PrismShape::PrismShape(unsigned int num_sides,
                       real_type    apothem,
                       real_type    rotate,
                       real_type    lo,
                       real_type    hi)
    : num_sides_(num_sides)
    , apothem_(apothem)
    , rotate_offset_(std::fmod(rotate, 1.0))
    , lo_(lo)
    , hi_(hi)
{
    CELER_VALIDATE(num_sides >= 3, << "Invalid number of sides " << num_sides);
    CELER_VALIDATE(lo <= hi,
                   << "Lower extent " << lo
                   << " must be less than upper extent " << hi << " for prism");
    CELER_VALIDATE(apothem > 0, << "Apothem for prism must be positive");
    CELER_VALIDATE(rotate >= 0 && rotate <= 1.0,
                   << "Rotation for prism must be in the range [0,1]");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string PrismShape::type() const
{
    return "prism";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool PrismShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief PrismShape volume
 */
real_type PrismShape::volume() const
{
    // See https://en.wikipedia.org/wiki/Regular_polygon#Area
    real_type xs_area = (num_sides_ * ipow<2>(apothem_))
                        * std::tan(pi / num_sides_);

    real_type infty = std::numeric_limits<real_type>::infinity();

    // Do not allow calculating volume for infinite prism shape
    CELER_VALIDATE(hi_ != infty && lo_ != -infty,
                   << "Volume calculation is not"
                      "allowed for infinite prism shape.");

    return xs_area * (hi_ - lo_);
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type PrismShape::inradius() const
{
    return apothem_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct surfaces for this shape
 */
void PrismShape::build(ShapeBuilder& build) const
{
    // Offset (if user offset is zero) is calculated to put a plane on the -y
    // face (sitting upright as visualized).. An offset of 1 produces a shape
    // congruent with an offset of zero, except that every face has an index
    // that's decremented by 1.
    const real_type offset
        = std::fmod(num_sides_ * 3.0 + 4 * rotate_offset_, 4.0) / 4.0;
    CELER_ASSERT(offset >= 0.0 && offset <= 1.0);

    // Fractions of a full rotation per side
    const real_type delta_turns = 1.0 / static_cast<real_type>(num_sides_);

    // Build prismatic sides
    for (auto n : range(num_sides_))
    {
        real_type       rot   = (n + offset) * delta_turns;
        const real_type theta = two_pi * rot;

        // Create a normal vector along the X axis, then rotate it through the
        // angle theta
        Real3 normal(0);
        normal[Axis::x] = std::cos(theta);
        normal[Axis::y] = std::sin(theta);

        // offset = normal * apotherm
        Real3 displacement(normal);
        displacement *= apothem_;

        build.plane(normal, displacement);
    }

    // Build top and bottom
    build.plane(Axis::z, pos, lo_);
    build.plane(Axis::z, neg, hi_);

    // Clip radial directions based on outer shape radius
    real_type radius = this->calc_circumradius();
    CELER_ASSERT(radius > apothem_);

    for (auto ax : {Axis::x, Axis::y})
    {
        build.clip_bbox_lower(ax, -radius);
        build.clip_bbox_upper(ax, radius);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Circumradius (outer radius) of the prism
 */
real_type PrismShape::calc_circumradius() const
{
    return apothem_ / std::cos(pi / num_sides_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
