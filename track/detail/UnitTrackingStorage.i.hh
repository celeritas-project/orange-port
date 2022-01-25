//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/UnitTrackingStorage.i.hh
 * \brief UnitTrackingStorage inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Macros.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Access all stored cells
 */
CELER_FORCEINLINE_FUNCTION const CellContainer&
                                 UnitTrackingStorage::cells() const
{
    return cells_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access all stored surfaces
 */
CELER_FORCEINLINE_FUNCTION const SurfaceContainer&
                                 UnitTrackingStorage::surfaces() const
{
    return surfaces_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find KD-tree cells at the given point
 */
auto UnitTrackingStorage::find_kdtree_cells(SpanConstReal3 pos) const
    -> SpanConstCellInt
{
    return make_span(kdtree_.find_volumes(pos));
}

//---------------------------------------------------------------------------//
/*!
 * Find cells connected to a given surface
 */
auto UnitTrackingStorage::get_connected_cells(SurfaceId id) const
    -> SpanConstCellInt
{
    CELER_EXPECT(id < connectivity_.size());
    return make_span(connectivity_[id.get()]);
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
