//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CuboidShape.cc
 * \brief CuboidShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CuboidShape.hh"

#include "base/GeometryUtils.hh"
#include "base/Range.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Construct from alternating lo/hi pairs
 */
CuboidShape CuboidShape::from_bounds(PairDbl x, PairDbl y, PairDbl z)
{
    return CuboidShape({x.first, y.first, z.first},
                       {x.second, y.second, z.second});
}

//---------------------------------------------------------------------------//
/*!
 * Construct a cuboid centered on the origin
 */
CuboidShape::CuboidShape(Real3 lo, Real3 hi)
    : lo_(std::move(lo)), hi_(std::move(hi))
{
    for (auto ax : range(def::END_XYZ))
    {
        CELER_VALIDATE(lo_[ax] <= hi_[ax],
                       << "Invalid cuboid " << to_cstring(ax) << " extents ["
                       << lo_[ax] << "," << hi_[ax] << "]");
    }
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Shape class name
 */
std::string CuboidShape::type() const
{
    return "cuboid";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool CuboidShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief CuboidShape volume
 */
real_type CuboidShape::volume() const
{
    return (hi_[0] - lo_[0]) * (hi_[1] - lo_[1]) * (hi_[2] - lo_[2]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Largest sphere radius that fits in this shape
 */
real_type CuboidShape::inradius() const
{
    return std::min(hi_[0] - lo_[0], std::min(hi_[1] - lo_[1], hi_[2] - lo_[2]))
           / 2;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct surfaces for this shape
 */
void CuboidShape::build(ShapeBuilder& build) const
{
    build.plane(Axis::x, pos, lo_[0]);
    build.plane(Axis::x, neg, hi_[0]);
    build.plane(Axis::y, pos, lo_[1]);
    build.plane(Axis::y, neg, hi_[1]);
    build.plane(Axis::z, pos, lo_[2]);
    build.plane(Axis::z, neg, hi_[2]);
}

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Make a cube shape
 *
 * \todo Possibly replace later with a 'cube' shape to make shape types better
 * reflect user input.
 */
std::shared_ptr<const Shape> make_cube(real_type lo, real_type hi)
{
    return std::make_shared<CuboidShape>(Real3(lo, lo, lo), Real3(hi, hi, hi));
}

//---------------------------------------------------------------------------//
} // namespace celeritas
