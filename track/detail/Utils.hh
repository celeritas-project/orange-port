//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/Utils.hh
 * \brief Utils class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "orange/surfaces/SurfaceContainer.hh"
#include "orange/Definitions.hh"
#include "base/OpaqueIdRange.hh"
#include "CellContainer.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
// Calculate senses inside a cell.
inline void calc_senses(const CellContainer::mapped_type& cell_def,
                        const SurfaceContainer&           surfaces,
                        SpanConstReal3                    pos,
                        SenseContainer*                   senses,
                        FaceId* result_face = nullptr);

//---------------------------------------------------------------------------//
// Calculate the surface normal
inline Real3 calc_normal(const SurfaceContainer& surfaces,
                         SpanConstReal3          pos,
                         SurfaceId               surface,
                         Sense                   sense);

//---------------------------------------------------------------------------//
// Find the face index of the given (optional) surface in the cell
inline FaceId
find_face(const span<const SurfaceId>& surfaces, SurfaceId surface);

//---------------------------------------------------------------------------//
// Loop over all faces in the given cell
inline OpaqueIdRange<FaceId>
all_faces(const CellContainer::mapped_type& cell_def);

//---------------------------------------------------------------------------//
// Calculate the bump distance for a position
inline real_type calc_bump(SpanConstReal3 pos);

//---------------------------------------------------------------------------//
//! Comparator for face/distance pairs
struct CloserFace
{
    bool operator()(const NextFace& a, const NextFace& b) const
    {
        return a.second > 0 && (b.second <= 0 || a.second < b.second);
    }
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
#include "Utils.i.hh"
//---------------------------------------------------------------------------//
