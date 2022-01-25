//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CylinderSegmentShape.cc
 * \brief CylinderSegmentShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CylinderSegmentShape.hh"

#include <cmath>

#include "base/Assert.hh"
#include "base/Constants.hh"
#include "base/GeometryUtils.hh"
#include "detail/ShapeBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * Constructor.
 */
CylinderSegmentShape::CylinderSegmentShape(real_type inner,
                                           real_type outer,
                                           real_type beg_rad,
                                           real_type angle,
                                           real_type lo,
                                           real_type hi)
    : inner_(inner)
    , outer_(outer)
    , lo_(lo)
    , hi_(hi)
    , beg_rad_(beg_rad)
    , end_rad_(beg_rad + angle)
{
    CELER_VALIDATE(lo <= hi,
                   << "Lower extent " << lo
                   << " must be less than upper extent " << hi
                   << " for cylinder segment");
    CELER_VALIDATE(outer > inner,
                   << "Outer radius " << outer
                   << " must be greater than inner radius " << inner
                   << " for cylinder segment");
    CELER_VALIDATE(inner >= 0,
                   << "Inner radius " << inner
                   << " must not be negative for cylinder segment");
    CELER_VALIDATE(angle <= constants::pi,
                   << "The angle subtended by the cylinder " << angle
                   << " must be less than pi (1/2 turn");
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Shape class name
 */
std::string CylinderSegmentShape::type() const
{
    return "cylinder_segment";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the shape is convex (no internal surface crossings)
 */
bool CylinderSegmentShape::is_convex() const
{
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Shape volume
 *
 * \verbatim
   V = (height) * (ring area) * (radial fraction)
     = (dz) * (pi (r_2^2 - r_1^2) * [(theta_2 - theta_1) / (2pi)]
 * \endverbatim
 */
real_type CylinderSegmentShape::volume() const
{
    return half * (hi_ - lo_) * (ipow<2>(outer_) - ipow<2>(inner_))
           * (end_rad_ - beg_rad_);
}

//---------------------------------------------------------------------------//
/*!
 * Largest sphere radius that fits in this shape
 */
real_type CylinderSegmentShape::inradius() const
{
    // Defined as the distance to maximum radii and half height
    // This assumes the maximum segment scenario of pi which involves the
    // origin
    real_type hh = half * (hi_ - lo_);
    return std::sqrt(ipow<2>(hh) + 2.0 * ipow<2>(outer_));
}

//---------------------------------------------------------------------------//
/*!
 * Construct surfaces for this shape
 */
void CylinderSegmentShape::build(ShapeBuilder& build) const
{
    // Build cylinders
    build.cyl(Axis::z, pos, inner_);
    build.cyl(Axis::z, neg, outer_);

    // Build bounding planes (above lower, below upper)
    build.plane(Axis::z, pos, lo_);
    build.plane(Axis::z, neg, hi_);

    // Build the azimuthal bounding surfaces for the cylindrical segment: The
    // normals are always defined *outward* because build.plane() for a
    // general plane defines the SENSE to be negative (in the opposite side of
    // the normal)
    Real3 origin(0.0, 0.0, 0.0);

    // Make outward normals using the azimuthal angles
    //   Begin: v = (cosphi, sinphi), n = (sinphi, -cosphi)
    //   End  : v = (cosphi, sinphi), n = (-sinphi, cosphi)
    Real3 n_beg(std::sin(beg_rad_), -std::cos(beg_rad_), 0.0);
    Real3 n_end(-std::sin(end_rad_), std::cos(end_rad_), 0.0);

    // Build the surfaces
    build.plane(n_beg, origin);
    build.plane(n_end, origin);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
