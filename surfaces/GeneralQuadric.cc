//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/GeneralQuadric.cc
 * \brief GeneralQuadric class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "GeneralQuadric.hh"

#include <cmath>
#include <iostream>
#include "Teuchos_BLAS.hpp"
#include "base/BLAS.hh"
#include "base/FixedViewArray.hh"
#include "base/VectorFunctions.hh"
#include "orange/Transform.hh"
#include "orange/BoundingBox.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * Construct with all coefficients
 */
GeneralQuadric::GeneralQuadric(const Real3& abc,
                               const Real3& def,
                               const Real3& ghi,
                               real_type    j)
    : a_(abc[0])
    , b_(abc[1])
    , c_(abc[2])
    , d_(def[0])
    , e_(def[1])
    , f_(def[2])
    , g_(ghi[0])
    , h_(ghi[1])
    , i_(ghi[2])
    , j_(j)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a new GeneralQuadric translated by some vector.
 *
 * See doc/quadric-transforms.pdf and nb/quadric-transform.ipynb
 */
GeneralQuadric GeneralQuadric::translated(const Transform& t) const
{
    CELER_EXPECT(!t.has_rotation());

    const real_type half_d = d_ / 2;
    const real_type half_e = e_ / 2;
    const real_type half_f = f_ / 2;

    // Linear components of the quadric matrix
    Real3 linear(g_ / 2, h_ / 2, i_ / 2);

    // Non-linear components of the quadric matrix
    SpaceMatrix nonl(Real3(a_, half_d, half_f),
                     Real3(half_d, b_, half_e),
                     Real3(half_f, half_e, c_));

    // Calculate q' = q - Q t
    Real3 newlinear(linear);
    blas().GEMV(Teuchos::NO_TRANS,
                3,
                3,
                -1.0,
                nonl[0].data(),
                3,
                t.translation().data(),
                1,
                1.0,
                newlinear.data(),
                1);

    // Create new quadric
    GeneralQuadric result(*this);

    // Update constant:
    // j' = j - t*(q' + q)
    result.j_ -= (t.translation()[0] * (linear[0] + newlinear[0])
                  + t.translation()[1] * (linear[1] + newlinear[1])
                  + t.translation()[2] * (linear[2] + newlinear[2]));

    // Update linear terms
    result.g_ = 2 * newlinear[0];
    result.h_ = 2 * newlinear[1];
    result.i_ = 2 * newlinear[2];
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Return a GeneralQuadric transformed by a rotate/translate
 *
 * See \c nb/quadric-transform.ipynb
 */
GeneralQuadric GeneralQuadric::transformed(const Transform& t) const
{
    CELER_EXPECT(t.has_rotation());

    using QuadVec = Array<real_type, 4>;
    using QuadMat = Array<QuadVec, 4>;

    // Build inverse transform matrix *COLUMN-MAJOR*
    QuadMat t_inv;
    {
        const Transform trans_inv = t.calc_inverse();
        t_inv[0][0]               = 1.0;
        for (int i = 1; i < 4; ++i)
        {
            t_inv[i][0] = 0.0;
        }

        const real_type* rot_iter = trans_inv.rotation_data();
        for (int j = 1; j < 4; ++j)
        {
            t_inv[0][j] = trans_inv.translation()[j - 1];
            for (int i = 1; i < 4; ++i)
            {
                t_inv[j][i] = *rot_iter++;
            }
        }
        CELER_ASSERT(rot_iter == trans_inv.rotation_data() + 9);
    }

    // Apply original quadric matrix to inverse transform
    QuadMat qprime;
    {
        QuadMat q(QuadVec(j_, g_ / 2, h_ / 2, i_ / 2),
                  QuadVec(0, a_, d_ / 2, f_ / 2),
                  QuadVec(0, 0, b_, e_ / 2),
                  QuadVec(0, 0, 0, c_));

        // Note: since the above input is column-major, upper triangular
        // becomes lower triangular
        blas().SYMM(Teuchos::LEFT_SIDE,
                    Teuchos::LOWER_TRI,
                    4,
                    4,
                    1.0,
                    q[0].data(),
                    4,
                    t_inv[0].data(),
                    4,
                    0,
                    qprime[0].data(),
                    4);
    }

    // Apply transpose of inverse transform to quadric
    QuadMat result;
    {
        blas().GEMM(Teuchos::TRANS,
                    Teuchos::NO_TRANS,
                    4,
                    4,
                    4,
                    1.0,
                    t_inv[0].data(),
                    4,
                    qprime[0].data(),
                    4,
                    0,
                    result[0].data(),
                    4);
    }

    // Create GQ from the result
    return GeneralQuadric(
        {result[1][1], result[2][2], result[3][3]},
        {2 * result[1][2], 2 * result[2][3], 2 * result[1][3]},
        {2 * result[0][1], 2 * result[0][2], 2 * result[0][3]},
        result[0][0]);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Describe to a stream.
 */
std::ostream& operator<<(std::ostream& os, const GeneralQuadric& s)
{
    os << "GQuadric: " << s.second() << ' ' << s.cross() << ' ' << s.first()
       << ' ' << s.zeroth();

    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
