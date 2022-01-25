//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/FoundCells.hh
 * \brief FoundCells class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Nemesis/containers/Span.hh"
#include "orange/Definitions.hh"
#include "CellContainer.hh"
#include "../Definitions.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Keep track of cell IDs during initialize for overlap checking.
 *
 * When finding a local cell, we need to ensure that only *one* cell
 * encloses the point in space. This class acts like a vector of cell IDs
 * where we only keep the two with the highest Z order. Insertion returns
 * 'true' with equal z order, and 'false' the first time a region with lower z
 * order is encountered.
 *
 * \note The z order checking relies on operating on cell IDs in descending
 * zorder. This is enforced by CellContainer.cc during construction. `insert`
 * has an extra assertion that the z order is only less or equal, never greater
 * (indicating unsortedness).
 */
class FoundCells
{
  public:
    //@{
    //! Public type aliases
    using size_type       = zorder_int;
    using SpanConstCellId = span<const VolumeId>;
    using SpanConstFaceId = span<const FaceId>;
    //@}

  public:
    // Constructor
    FoundCells() = default;

    // Add a new cell; return success (equal z order)
    inline bool insert(VolumeId id, zorder_int zorder, FoundFace face);

    //! Whether no cells were found
    bool empty() const { return size_ == size_type(0); }

    //! Number of found cells
    size_type size() const { return size_; }

    //! Z order of found cells (NONE if no found cells)
    zorder_int zorder() const { return cur_zorder_; }

    // View to the highest-precedence cells
    inline SpanConstCellId cell_view() const;

    // View to the face corresponding to found cell
    inline FaceId face() const;

    // Sense of the found cell
    inline Sense sense() const;

  private:
    //// DATA ////
    static constexpr int max_size = 2;

    zorder_int cur_zorder_ = zorder_int(ZOrder::invalid);
    size_type  size_       = 0;
    VolumeId   cell_storage_[max_size];
    FaceId     face_;
    Sense      sense_{};
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "FoundCells.i.hh"
//---------------------------------------------------------------------------//
