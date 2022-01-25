//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/Universe.cc
 * \brief Universe class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "Universe.hh"

#include <algorithm>
#include "base/Assert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with parameters.
 *
 * All specified daughters must be assigned and have a cell ID within the
 * allowed range.
 */
Universe::Universe(Params params)
    : id_(std::move(params.id))
    , tracker_(std::move(params.tracker))
    , daughters_(std::move(params.daughters))
{
    CELER_EXPECT(id_);
    CELER_EXPECT(tracker_);
    CELER_EXPECT(tracker_->num_volumes() == daughters_.size()
                 || daughters_.empty());
}

//---------------------------------------------------------------------------//
} // namespace celeritas
