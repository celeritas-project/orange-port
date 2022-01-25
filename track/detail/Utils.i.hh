//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/Utils.i.hh
 * \brief Utils inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <algorithm>
#include <cmath>
#include "orange/surfaces/SurfaceAction.hh"
#include "orange/Fuzziness.hh"
#include "SurfaceFunctors.hh"
#include "LogicEvaluator.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Calculate senses inside a cell.
 */
void calc_senses(const CellContainer::mapped_type& cell_def,
                 const SurfaceContainer&           surfaces,
                 SpanConstReal3                    pos,
                 SenseContainer*                   senses,
                 FaceId*                           result_face)
{
    CELER_EXPECT(senses);

    // Build a functor to calculate the sense of a surface ID given the current
    // state position
    auto calc_sense = make_surface_action(surfaces, detail::CalcSense{pos});

    // Resulting "on" face, may already be set
    FaceId on_face = result_face ? *result_face : FaceId{};

    // Fill the temp logic vector with values for all surfaces in the cell
    senses->resize(cell_def.faces.size());
    auto first      = cell_def.faces.begin();
    auto last       = cell_def.faces.end();
    auto sense_iter = senses->begin();
    for (auto surf_iter = first; surf_iter != last; ++surf_iter, ++sense_iter)
    {
        SignedSense ss = calc_sense(*surf_iter);
        *sense_iter    = to_sense(ss); // On or outside -> outside
        if (!on_face && ss == SignedSense::on)
        {
            // First face that we're exactly on
            on_face = FaceId(surf_iter - first);
        }
    }

    if (result_face)
    {
        // Save result face
        *result_face = on_face;
        CELER_ENSURE(!*result_face || *result_face < cell_def.faces.size());
    }
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the surface normal.
 */
Real3 calc_normal(const SurfaceContainer& surfaces,
                  SpanConstReal3          pos,
                  SurfaceId               surface,
                  Sense                   sense)
{
    CELER_EXPECT(surface);
    auto calc_impl = make_surface_action(surfaces, CalcNormal{pos});

    Real3 normal = calc_impl(surface);
    // TODO: re-enable and update tests
#if 0
    if (sense == Sense::neg)
    {
        // Flip sense so it's "outward" just for jollies
        normal *= -1;
    }
#endif

    return normal;
}

//---------------------------------------------------------------------------//
/*!
 * Find the index of a surface inside a cell.
 *
 * The surface ID stored on the state should be local to the cell; the
 * resulting *cell* surface ID is relative to the surfaces for that cell.
 */
FaceId find_face(const span<const SurfaceId>& surfaces, SurfaceId surface)
{
    FaceId result;
    if (surface)
    {
        // Find index of surface in sorted list
        auto iter = std::lower_bound(surfaces.begin(), surfaces.end(), surface);
        if (iter != surfaces.end() && *iter == surface)
        {
            result = FaceId(iter - surfaces.begin());
        }
    }
    CELER_ENSURE(!result || result < surfaces.size());
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Loop over all faces in the given cell
 */
OpaqueIdRange<FaceId> all_faces(const CellContainer::mapped_type& cell_def)
{
    return id_range(FaceId(cell_def.faces.size()));
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the bump distance for a position
 */
real_type calc_bump(SpanConstReal3 pos)
{
    const auto& fuzz = fuzziness();

    // Find the greatest of "absolute" bump with the largest-magnitude relative
    // bumps
    real_type dist = fuzz.bump_abs();
    for (real_type p : pos)
    {
        dist = std::max(dist, fuzz.bump_rel() * std::fabs(p));
    }
    CELER_ENSURE(dist > 0);
    return dist;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
