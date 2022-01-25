//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/Transform.cc
 * \brief Transform class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Transform.hh"

#include <array>
#include "Teuchos_BLAS.hpp"
#include "base/Assert.hh"
#include "base/Definitions.hh"
#include "base/BLAS.hh"
#include "orange/Transform.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
void transform_to_daughter(SpanConstReal3   parent_pos,
                           SpanConstReal3   parent_dir,
                           const Transform& transform,
                           SpanReal3        pos,
                           SpanReal3        dir)
{
    std::array<std::array<real_type, 3>, 2> pos_dir;

    // Translate before applying rotation
    const auto& translation = transform.translation();
    for (auto i : {Axis::x, Axis::y, Axis::z})
    {
        pos_dir[0][i] = parent_pos[i] - translation[i];
        pos_dir[1][i] = parent_dir[i];
    }

    // Only apply rotation if there is a non-identity matrix
    if (transform.has_rotation())
    {
        // Conveniently, BLAS wants things in COLUMN_MAJOR order, so we can
        // simply concatenate our vectors to construct the matrix

        constexpr real_type alpha    = 1.0;
        constexpr int       rot_cols = 3;
        const real_type*    rot      = transform.rotation_data();

        constexpr real_type                     beta = 0.0;
        std::array<std::array<real_type, 3>, 2> dst;

        // Use BLAS to do the matrix-matrix multiply
        constexpr int stride = 3;
        blas().GEMM(Teuchos::TRANS,
                    Teuchos::NO_TRANS,
                    rot_cols,
                    dst.size(),
                    dst[0].size(),
                    alpha,
                    rot,
                    stride,
                    pos_dir[0].data(),
                    stride,
                    beta,
                    dst[0].data(),
                    stride);
        std::copy(
            dst[0].data(), dst[1].data() + dst[1].size(), pos_dir[0].data());
    }

    // Copy the results back
    std::copy(pos_dir[0].begin(), pos_dir[0].end(), pos.data());
    std::copy(pos_dir[1].begin(), pos_dir[1].end(), dir.data());
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
