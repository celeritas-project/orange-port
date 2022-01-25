//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/ObjectMetadata.i.hh
 * \brief ObjectMetadata inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Provenance (originating file/line; always set)
 */
auto ObjectMetadata::provenance() const -> const SPConstProvenance&
{
    CELER_EXPECT(*this);
    return data_.provenance;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
