//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ParallelepipedShape.cc
 * \brief ParallelepipedShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ParallelepipedShape.hh"

#include <cmath>
#include <string>

#include "base/Constants.hh"
#include "base/Assert.hh"
#include "base/Future.hh"
#include "base/Range.hh"
#include "base/VectorFunctions.hh"
#include "orange/surfaces/GeneralQuadric.hh"
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

using Axis::x;
using Axis::y;
using Axis::z;
using def::I;
using def::J;
using def::K;

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Parallelepiped Constructor.
 */
ParallelepipedShape::ParallelepipedShape(Real3     length,
                                         real_type psi,
                                         real_type theta,
                                         real_type phi)
{
    using soft_equiv;
    using vector_magnitude;
    using constants::pi;
    constexpr real_type half_pi = pi / 2.0;

    CELER_VALIDATE(length[X] > 0.0,
                   << "ppiped x length " << length[X] << " must be positive");
    CELER_VALIDATE(length[Y] > 0.0,
                   << "ppiped y length " << length[Y] << " must be positive");
    CELER_VALIDATE(length[Z] > 0.0,
                   << "ppiped z length " << length[Z] << " must be positive");

    // Must be less than half_pi degrees (positive volume), i.e. pi/2
    CELER_VALIDATE(psi >= 0.0 && psi < half_pi,
                   << "angle between x face and y-axis " << psi
                   << " must be [0,half_pi)");
    CELER_VALIDATE(theta >= 0.0 && theta < half_pi,
                   << "angle between y face and z-axis " << theta
                   << " must be [0,half_pi)");
    CELER_VALIDATE(phi >= 0.0 && phi < half_pi,
                   << "angle between x face and x-axis " << phi
                   << " must be [0,half_pi)");

    basis_[I][X] = length[X];

    basis_[J][X] = std::sin(psi) * length[Y];
    basis_[J][Y] = std::cos(psi) * length[Y];

    const real_type stheta = std::sin(theta);
    const real_type ctheta = std::cos(theta);

    basis_[K][X] = stheta * std::cos(phi) * length[Z];
    basis_[K][Y] = stheta * std::sin(phi) * length[Z];
    basis_[K][Z] = ctheta * length[Z];

    CELER_ENSURE(soft_equiv(vector_magnitude(basis_[I]), length[X]));
    CELER_ENSURE(soft_equiv(vector_magnitude(basis_[J]), length[Y]));
    CELER_ENSURE(soft_equiv(vector_magnitude(basis_[K]), length[Z]));
}
//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string ParallelepipedShape::type() const
{
    return "parallelepiped";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool ParallelepipedShape::is_convex() const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief ParallelepipedShape volume
 */
real_type ParallelepipedShape::volume() const
{
    // Triple product of scaled basis vectors
    real_type vol = dot_product(basis_[K], cross_product(basis_[I], basis_[J]));

    CELER_ENSURE(vol > 0.0);
    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type ParallelepipedShape::inradius() const
{
    // Unknown
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void ParallelepipedShape::build(ShapeBuilder& build) const
{
    // Create the x plane pair
    Real3 xp_normal = cross_product(basis_[J], basis_[K]);
    build.plane(-xp_normal, Real3(0, 0, 0));
    build.plane(xp_normal, basis_[I]);

    // Create the y plane pair
    Real3 yp_normal = cross_product(basis_[K], basis_[I]);
    build.plane(-yp_normal, Real3(0, 0, 0));
    build.plane(yp_normal, basis_[J]);

    Real3 upper_vertex = calc_upper_vertex();

    // bottom and top z planes
    build.plane(Z, pos, 0);
    build.plane(Z, neg, upper_vertex[Z]);

    // Set lower side of bounding box
    build.clip_bbox_lower(X, 0);
    build.clip_bbox_lower(Y, 0);
    build.clip_bbox_lower(Z, 0);

    // Set upper side of bounding box
    build.clip_bbox_upper(X, upper_vertex[X]);
    build.clip_bbox_upper(Y, upper_vertex[Y]);
    build.clip_bbox_upper(Z, upper_vertex[Z]);
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
/*!
 * Calculate the upper vertex of the ppiped.
 *
 * This is just the vector sum of the three scaled basis functions.
 */
Real3 ParallelepipedShape::calc_upper_vertex() const
{
    Real3 upper_vertex = basis_[I];
    upper_vertex += basis_[J];
    upper_vertex += basis_[K];
    return upper_vertex;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
