//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/UnitTrackingStorage.hh
 * \brief UnitTrackingStorage class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Nemesis/containers/Span.hh"
#include "base/FlatTable.hh"
#include "orange/kdtree/KDTree.hh"
#include "orange/surfaces/SurfaceContainer.hh"
#include "CellContainer.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Container for surface/cell data needed for unit trackers.
 */
class UnitTrackingStorage
{
  public:
    //@{
    //! Public type aliases
    using volume_int       = geometria::cell_type;
    using SpanConstCellInt = span<const volume_int>;
    //@}

  public:
    // Construct from surfaces and regions
    UnitTrackingStorage(SurfaceContainer               surfaces,
                        const std::vector<UnitRegion>& regions);

    //// SEARCH ////

    // Find KD-tree cells at the given point
    inline SpanConstCellInt find_kdtree_cells(SpanConstReal3 pos) const;

    // Find cells connected to a given surface
    inline SpanConstCellInt get_connected_cells(SurfaceId id) const;

    //// ACCESSORS ////

    //! Number of cells
    size_type num_volumes() const { return cells_.size(); }

    //! Number of surfaces
    size_type num_surfaces() const { return surfaces_.size(); }

    //! Access all cells
    inline const CellContainer& cells() const;

    //! Access all surfaces
    inline const SurfaceContainer& surfaces() const;

  private:
    //// TYPEDEFS ////

    using TableCellId = FlatTable<volume_int>;

    //// DATA ////

    //! Compressed vector of surface definitions
    SurfaceContainer surfaces_;

    //! Compressed vector of flagged cell definitions
    detail::CellContainer cells_;

    //! KD-tree acceleration structure
    geometria::KDTree kdtree_;

    //! Surface ID -> [Cell ID, ...] connectivity after build
    TableCellId connectivity_;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "UnitTrackingStorage.i.hh"
//---------------------------------------------------------------------------//
