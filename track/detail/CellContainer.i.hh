//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/CellContainer.i.hh
 * \brief CellContainer inline method definitions
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
 * Get a cell reference
 */
CELER_FORCEINLINE_FUNCTION auto CellContainer::operator[](VolumeId id) const
    -> const mapped_type&
{
    CELER_EXPECT(id < this->size());
    return cells_[id.get()];
}

//---------------------------------------------------------------------------//
/*!
 * Set flags for a cell
 */
void CellContainer::set_flags(VolumeId id, flag_int flags)
{
    CELER_EXPECT(id < this->size());
    CELER_EXPECT(flags >= 0);
    cells_[id.get()].flags = flags;
}

//---------------------------------------------------------------------------//
/*!
 * Set number of intersections for a cell
 */
void CellContainer::set_num_intersections(VolumeId id, logic_int val)
{
    CELER_EXPECT(id < this->size());
    CELER_EXPECT(val >= 0);
    cells_[id.get()].num_intersections = val;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
