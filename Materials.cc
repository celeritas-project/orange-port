//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file Materials.cc
 * \brief Materials class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Materials.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with matids and indexer.
 */
Materials::Materials(VecMatid matids, SPConstIndexer indexer)
    : matids_(std::move(matids)), indexer_(std::move(indexer))
{
    CELER_EXPECT(indexer_);
    CELER_EXPECT(matids_.size() == indexer_->num_volumes());
}

//---------------------------------------------------------------------------//
} // namespace celeritas
