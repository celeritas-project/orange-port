//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/SimpleQuadric.cc
 * \brief SimpleQuadric class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SimpleQuadric.hh"

#include <cmath>
#include <iostream>
#include "base/FixedViewArray.hh"
#include "base/VectorFunctions.hh"
#include "orange/Transform.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "GeneralQuadric.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * Whether the given quadric can be simplified to this shape
 */
bool SimpleQuadric::can_simplify(const GeneralQuadric& gq)
{
    const real_type tol   = fuzziness().surface_simplification_abs();
    auto            cross = gq.cross();
    return std::fabs(cross[0]) < tol && std::fabs(cross[1]) < tol
           && std::fabs(cross[2]) < tol;
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * Construct at the origin
 */
SimpleQuadric::SimpleQuadric(const Real3& abc, const Real3& def, real_type g)
    : a_(abc[0])
    , b_(abc[1])
    , c_(abc[2])
    , d_(def[0])
    , e_(def[1])
    , f_(def[2])
    , g_(g)
{
    // Ill-defined if all non-constants are zero
    CELER_EXPECT(a_ != 0 || b_ != 0 || c_ != 0 || d_ != 0 || e_ != 0
                 || f_ != 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with all coefficients
 */
SimpleQuadric::SimpleQuadric(const Real3& abc,
                             const Real3& def,
                             real_type    g,
                             const Real3& origin)
    : SimpleQuadric(abc, def, g)
{
    const Real3 orig_def = make_vector(def);

    // Expand out origin into the other terms
    d_ -= 2 * a_ * origin[0];
    e_ -= 2 * b_ * origin[1];
    f_ -= 2 * c_ * origin[2];

    g_ += a_ * origin[0] * origin[0];
    g_ += b_ * origin[1] * origin[1];
    g_ += c_ * origin[2] * origin[2];

    g_ -= 2 * orig_def[0] * origin[0];
    g_ -= 2 * orig_def[1] * origin[1];
    g_ -= 2 * orig_def[2] * origin[2];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from a degenerate GQ
 *
 * This just deletes the cross terms.
 */
SimpleQuadric::SimpleQuadric(const GeneralQuadric& gq)
    : SimpleQuadric(gq.second(), gq.first(), gq.zeroth())
{
    CELER_EXPECT(can_simplify(gq));
}

//---------------------------------------------------------------------------//
/*!
 * Return a new SimpleQuadric translated by some vector
 */
SimpleQuadric SimpleQuadric::translated(const Transform& t) const
{
    CELER_EXPECT(!t.has_rotation());

    return SimpleQuadric(
        this->second(), this->first(), this->zeroth(), t.translation());
}

//---------------------------------------------------------------------------//
/*!
 * Return a SimpleQuadric transformed by a rotate/translate
 */
GeneralQuadric SimpleQuadric::transformed(const Transform& t) const
{
    CELER_EXPECT(t.has_rotation());

    return GeneralQuadric(
               this->second(), {0, 0, 0}, this->first(), this->zeroth())
        .transformed(t);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Describe to a stream.
 */
std::ostream& operator<<(std::ostream& os, const SimpleQuadric& s)
{
    os << "SQuadric: " << s.second() << ' ' << s.first() << ' ' << s.zeroth();
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
