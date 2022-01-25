//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CylAligned.t.hh
 * \brief CylAligned template method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "CylAligned.hh"

#include <iostream>
#include <cmath>
#include "base/FixedViewArray.hh"
#include "orange/Transform.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "SimpleQuadric.hh"
#include "GeneralQuadric.hh"
#include "CylCentered.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * Whether the given quadric can be simplified to this shape.
 *
 * A cylindrical quadric has no cross terms, a zero second-order term along
 * this axis, and equal off-axis second-order terms.
 *
 * This implementation requires that the second-order terms be scaled to 1.
 *
 * Note that this will *not* allow simplification of an "inverse" cylinder
 * (an SQ defined to include the half-space outside the cylinder boundary). You
 * will first need to invert the coefficients (and the sense if used in that
 * context).
 */
template<Axis T>
bool CylAligned<T>::can_simplify(const SimpleQuadric& sq)
{
    auto            second = sq.second();
    const real_type tol    = fuzziness().surface_simplification_abs();
    return std::fabs(second[t_index()]) < tol
           && std::fabs(second[u_index()] - 1.0) < tol
           && std::fabs(second[v_index()] - 1.0) < tol;
}

//---------------------------------------------------------------------------//
// METHODS
//---------------------------------------------------------------------------//
/*!
 * User-accessible construction
 */
template<Axis T>
CylAligned<T>::CylAligned(Real3 origin, real_type radius)
    : radius_sq_ipow<2>(radius)
{
    CELER_EXPECT(radius >= 0);

    // Set off-axis coordinates
    origin_u_ = origin[u_index()];
    origin_v_ = origin[v_index()];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from a degenerate quadric
 */
template<Axis T>
CylAligned<T>::CylAligned(const SimpleQuadric& sq)
{
    CELER_EXPECT(can_simplify(sq));

    // Calc origin from first-order coefficients
    auto first = sq.first();
    origin_u_  = first[u_index()] * (-half);
    origin_v_  = first[v_index()] * (-half);
    radius_sq_ = -sq.zeroth() + ipow<2>(origin_u_) + ipow<2>(origin_v_);
    CELER_ENSURE(radius_sq_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * Construct from a centered ortho cyl
 */
template<Axis T>
CylAligned<T>::CylAligned(const CylCentered<T>& coc)
    : origin_u_(0), origin_v_(0), radius_sq_(coc.radius_sq())
{
    CELER_ENSURE(radius_sq_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * Return a new cylinder translated by some vector
 */
template<Axis T>
CylAligned<T> CylAligned<T>::translated(const Transform& t) const
{
    CELER_EXPECT(!t.has_rotation());

    // Create a copy
    CylAligned<T> result(*this);

    // Set off-axis coordinates
    result.origin_u_ += t.translation()[u_index()];
    result.origin_v_ += t.translation()[v_index()];

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Return a cylinder transformed by a rotate/translate
 *
 * We convert this to a general quadric, which is then tranformed.. The GQ
 * coefficients are found by expanding the cylinder's coefficients. For the
 * case of a cylinder parallel to the axis, this is \f[
  y^2 + z^2 - 2y_0 y - 2z_0 z + (y_0^2 + z_0^2 - R^2) = 0
  \f]
 */
template<Axis T>
GeneralQuadric CylAligned<T>::transformed(const Transform& t) const
{
    CELER_EXPECT(t.has_rotation());

    // First, second, and zeroth coefficients
    Real3     second;
    Real3     first;
    real_type zeroth;

    second[t_index()] = 0;
    second[u_index()] = 1;
    second[v_index()] = 1;

    first[t_index()] = 0;
    first[u_index()] = -2 * origin_u_;
    first[v_index()] = -2 * origin_v_;

    zeroth = -radius_sq_ + ipow<2>(origin_u_) + ipow<2>(origin_v_);

    return GeneralQuadric(make_fixed_view(second),
                          {0, 0, 0},
                          make_fixed_view(first),
                          zeroth)
        .transformed(t);
}

//---------------------------------------------------------------------------//
/*!
 * Clip a bounding box to this surface
 */
template<Axis T>
void CylAligned<T>::clip(Sense sense, BoundingBox& bbox) const
{
    if (sense == neg)
    {
        real_type radius = std::sqrt(radius_sq_);
        // Increase very slightly (bbox_rel is BVH fuzziness ~ 1e-7; we only
        // need to avoid clipping the actual cylinder edges..)
        radius *= 1 + fuzziness().bbox_rel() * fuzziness().bbox_rel();

        bbox.clip_lower(u_index(), origin_u_ - radius);
        bbox.clip_upper(u_index(), origin_u_ + radius);
        bbox.clip_lower(v_index(), origin_v_ - radius);
        bbox.clip_upper(v_index(), origin_v_ + radius);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Output a very short description of this sufrace
 */
template<Axis T>
std::ostream& operator<<(std::ostream& os, const CylAligned<T>& s)
{
    os << "Cyl " << to_cstring(T) << ": r=" << std::sqrt(s.radius_sq()) << ", "
       << to_cstring(CylAligned<T>::U) << '=' << s.origin_u() << ", "
       << to_cstring(CylAligned<T>::V) << '=' << s.origin_v();
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
