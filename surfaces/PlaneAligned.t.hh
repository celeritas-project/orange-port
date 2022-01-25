//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/PlaneAligned.t.hh
 * \brief PlaneAligned template method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "PlaneAligned.hh"

#include <iostream>
#include "orange/Transform.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "Plane.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * Whether the given quadric can be simplified to this shape.
 *
 * Note that for an ortho plane, the normal *must* be in the positive
 * half-space! Logic in the surface insertion utility should take care of this
 * flipping.
 */
template<Axis T>
bool PlaneAligned<T>::can_simplify(const Plane& p)
{
    const real_type tol = fuzziness().surface_simplification_abs();
    return std::fabs(p.normal()[t_index()] - 1.0) < tol;
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * User-accessible construction from position
 */
template<Axis T>
PlaneAligned<T>::PlaneAligned(real_type position) : position_(position)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from a degenerate plane
 */
template<Axis T>
PlaneAligned<T>::PlaneAligned(const Plane& p)
{
    CELER_EXPECT(can_simplify(p));

    // Change from "displacement" to "intersection", adjusting the location
    // slightly in case axis_normal isn't exactly 1.
    const real_type orig_norm = p.normal()[t_index()];
    position_                 = p.displacement() / orig_norm;
}

//---------------------------------------------------------------------------//
/*!
 * Return a new plane translated by some vector
 */
template<Axis T>
PlaneAligned<T> PlaneAligned<T>::translated(const Transform& t) const
{
    CELER_EXPECT(!t.has_rotation());
    real_type normal_translation = t.translation()[t_index()];
    return PlaneAligned<T>(position_ + normal_translation);
}

//---------------------------------------------------------------------------//
/*!
 * Return a plane transformed by a rotate/translate
 */
template<Axis T>
Plane PlaneAligned<T>::transformed(const Transform& t) const
{
    CELER_EXPECT(t.has_rotation());

    // Get the normal for this plane
    Real3 rn;
    rn[t_index()] = 1.0;

    // Get the original point on the plane
    Real3 rp(0, 0, 0);
    rp[t_index()] = position_;

    // Rotate the normal and transform the point
    t.rotate_to_parent(rn);
    t.daughter_to_parent(rp);

    // Build the general plane corresponding to this transform
    return Plane(rn, dot_product(rn, rp));
}

//---------------------------------------------------------------------------//
/*!
 * Clip a bounding box to this surface
 */
template<Axis T>
void PlaneAligned<T>::clip(Sense sense, BoundingBox& bbox) const
{
    if (sense == neg)
    {
        bbox.clip_upper(t_index(), position_);
    }
    else
    {
        bbox.clip_lower(t_index(), position_);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Output a very short description of this sufrace
 */
template<Axis T>
std::ostream& operator<<(std::ostream& os, const PlaneAligned<T>& s)
{
    os << "Plane: " << to_cstring(T) << '=' << s.position();
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
