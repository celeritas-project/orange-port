//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/FoundCells.i.hh
 * \brief FoundCells inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Assert.hh"
#include "base/Macros.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Add a new cell; return success (equal z order but not overflow)
 *
 * When max_size equal-zorder cells have been found, we return false.
 *
 * It's OK to only have a single face and sense, since having multiple valid
 * cells at the same Z order is an error
 */
CELER_FORCEINLINE_FUNCTION bool
FoundCells::insert(VolumeId id, zorder_int zorder, FoundFace found)
{
    CELER_EXPECT(id);
    CELER_EXPECT(zorder <= cur_zorder_ || this->empty());
    CELER_EXPECT(found);

    if (zorder < cur_zorder_)
    {
        // New cell is hidden by higher-priority cells
        return false;
    }
    else if (NEMESIS_UNLIKELY(this->size() >= this->max_size))
    {
        // Capacity is exceeded
        return false;
    }

    // Save zorder
    cur_zorder_ = zorder;

    // Save cell, face, sense
    cell_storage_[size_] = id;
    face_                = found.face;
    sense_               = found.sense;

    // Update size
    ++size_;

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * View to the highest-precedence cells
 */
CELER_FORCEINLINE_FUNCTION auto FoundCells::cell_view() const -> SpanConstCellId
{
    return {cell_storage_, size_};
}

//---------------------------------------------------------------------------//
/*!
 * Face corresponding to found cell.
 */
FaceId FoundCells::face() const
{
    CELER_EXPECT(this->size() == 1);
    return face_;
}

//---------------------------------------------------------------------------//
/*!
 * Sense corresponding to found face, if present (unspecified if not).
 */
Sense FoundCells::sense() const
{
    CELER_EXPECT(this->size() == 1);
    return sense_;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
