//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/ConeAligned.t.hh
 * \brief ConeAligned template method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "ConeAligned.hh"

#include <iostream>
#include <cmath>
#include "base/FixedViewArray.hh"
#include "orange/Transform.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "GeneralQuadric.hh"
#include "SimpleQuadric.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * Whether the given quadric can be simplified to this shape
 *
 * A cone has a negative second-moment along this axis, no cross terms, and
 * equal non-axis second moments.
 *
 * This implementation requires that the second-order terms be scaled to 1.
 */
template<Axis T>
bool ConeAligned<T>::can_simplify(const SimpleQuadric& gq)
{
    const real_type tol    = fuzziness().surface_simplification_abs();
    auto            second = gq.second();
    return second[t_index()] < 0 && std::fabs(second[u_index()] - 1.0) < tol
           && std::fabs(second[v_index()] - 1.0) < tol;
}

//---------------------------------------------------------------------------//
// METHODS
//---------------------------------------------------------------------------//
/*!
 * User-accessible construction
 */
template<Axis T>
ConeAligned<T>::ConeAligned(const Real3& origin, real_type tangent)
    : origin_(make_vector(origin)), tsq_ipow<2>(tangent)
{
    CELER_VALIDATE(tangent >= 0,
                   << "Cone tangent must be nonnegative, but it is "
                   << tangent);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from a degenerate SQ
 */
template<Axis T>
ConeAligned<T>::ConeAligned(const SimpleQuadric& sq)
{
    CELER_EXPECT(can_simplify(sq));

    // Calculate t^2 from the second-order term
    tsq_ = -sq.second()[t_index()];

    // Calc origin from first-order coefficients
    origin_[t_index()] = sq.first()[t_index()] / (2 * tsq_);
    origin_[u_index()] = sq.first()[u_index()] * (-half);
    origin_[v_index()] = sq.first()[v_index()] * (-half);

    CELER_ENSURE(soft_equiv(-tsq_ * origin_[t_index()] * origin_[t_index()]
                                + (origin_[u_index()] * origin_[u_index()])
                                + (origin_[v_index()] * origin_[v_index()]),
                            sq.zeroth(),
                            1.e-6));
}

//---------------------------------------------------------------------------//
/*!
 * Return a new cone translated by some vector
 */
template<Axis T>
ConeAligned<T> ConeAligned<T>::translated(const Transform& t) const
{
    CELER_EXPECT(!t.has_rotation());

    // Create a copy
    ConeAligned<T> result(*this);

    // Update origin
    result.origin_ += t.translation();

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Return a cone transformed by a rotate/translate
 *
 * A cone parallel to the X axis is represented as a general quadric:
 * \f[
    (-t^2) x^2 + y^2 + z^2 + (2t^2 x_0) x + (-2y_0) y + (-2z_0) z
    + (-t^2 x_0^2 + y_0^2 + z_0^2) = 0
   \f]
 */
template<Axis T>
GeneralQuadric ConeAligned<T>::transformed(const Transform& t) const
{
    CELER_EXPECT(t.has_rotation());

    // First, second, and zeroth coefficients
    Real3     second;
    Real3     first;
    real_type zeroth;

    second[t_index()] = -tsq_;
    second[u_index()] = 1;
    second[v_index()] = 1;

    first[t_index()] = 2 * tsq_ * origin_[t_index()];
    first[u_index()] = -2 * origin_[u_index()];
    first[v_index()] = -2 * origin_[v_index()];

    zeroth = -tsq_ * origin_[t_index()] * origin_[t_index()]
             + (origin_[u_index()] * origin_[u_index()])
             + (origin_[v_index()] * origin_[v_index()]);

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
void ConeAligned<T>::clip(Sense, BoundingBox&) const
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output a very short description of this sufrace
 */
template<Axis T>
std::ostream& operator<<(std::ostream& os, const ConeAligned<T>& s)
{
    os << "Cone " << to_cstring(T) << ": tangent=" << std::sqrt(s.tangent_sq())
       << " at " << s.origin();
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
