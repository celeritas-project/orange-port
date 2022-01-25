//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CylCentered.t.hh
 * \brief CylCentered template method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "CylCentered.hh"

#include <iostream>
#include <cmath>
#include "orange/Transform.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "CylAligned.hh"
#include "GeneralQuadric.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * Whether the given ortho cyl can be simplified to this shape.
 *
 * This uses the magnitude of the radius to determine whether the origin is
 * soft equivalent to zero.
 */
template<Axis T>
bool CylCentered<T>::can_simplify(const CylAligned<T>& oc)
{
    const real_type tol = fuzziness().surface_simplification_abs()
                          * oc.radius_sq();
    return (oc.origin_u() * oc.origin_u() < tol)
           && (oc.origin_v() * oc.origin_v() < tol);
}

//---------------------------------------------------------------------------//
// METHODS
//---------------------------------------------------------------------------//
/*!
 * User-accessible construction
 */
template<Axis T>
CylCentered<T>::CylCentered(real_type radius) : radius_sq_ipow<2>(radius)
{
    CELER_EXPECT(radius >= 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from a degenerate GQ
 */
template<Axis T>
CylCentered<T>::CylCentered(const CylAligned<T>& oc)
{
    CELER_EXPECT(can_simplify(oc));

    // Calc origin from first-order coefficients
    radius_sq_ = oc.radius_sq();
    CELER_ENSURE(radius_sq_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * Return a new cylinder translated by some vector
 */
template<Axis T>
CylAligned<T> CylCentered<T>::translated(const Transform& t) const
{
    CELER_EXPECT(!t.has_rotation());

    // Return a translated ortho cyl
    return CylAligned<T>(*this).translated(t);
}

//---------------------------------------------------------------------------//
/*!
 * Return a cylinder transformed by a rotate/translate
 *
 * This relies on the ortho cyl
 */
template<Axis T>
GeneralQuadric CylCentered<T>::transformed(const Transform& t) const
{
    CELER_EXPECT(t.has_rotation());

    // Return a transformed ortho cyl
    return CylAligned<T>(*this).transformed(t);
}

//---------------------------------------------------------------------------//
/*!
 * Clip a bounding box to this surface
 */
template<Axis T>
void CylCentered<T>::clip(Sense sense, BoundingBox& bbox) const
{
    if (sense == neg)
    {
        real_type radius = std::sqrt(radius_sq_);
        // Increase very slightly (bbox_rel is BVH fuzziness ~ 1e-7; we only
        // need to avoid clipping the actual cylinder edges..)
        radius *= 1 + fuzziness().bbox_rel() * fuzziness().bbox_rel();

        bbox.clip_lower(u_index(), -radius);
        bbox.clip_upper(u_index(), radius);
        bbox.clip_lower(v_index(), -radius);
        bbox.clip_upper(v_index(), radius);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Output a very short description of this sufrace
 */
template<Axis T>
std::ostream& operator<<(std::ostream& os, const CylCentered<T>& s)
{
    os << "Cyl " << to_cstring(T) << ": r=" << std::sqrt(s.radius_sq());
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
