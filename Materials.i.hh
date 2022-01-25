//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file Materials.i.hh
 * \brief Materials inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "orange/track/StateView.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Get matid from a cell ID.
 *
 * Needed for current Shift geometry interface
 */
auto Materials::matid(cell_type volid) const -> matid_type
{
    CELER_EXPECT(volid < matids_.size());
    return matids_[volid];
}

//---------------------------------------------------------------------------//
/*!
 * Get the material ID from the current state.
 *
 * The ORANGEGeometry interface class can cache this as necessary.
 */
auto Materials::matid(const ORANGEState& state) const -> matid_type
{
    const StateView state_view(state.tracking);

    auto volid
        = indexer_->global_cell(state_view.universe(), state_view.volume());
    return this->matid(volid);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
