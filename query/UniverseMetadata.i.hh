//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/UniverseMetadata.i.hh
 * \brief UniverseMetadata inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <utility>
#include "detail/DescribeMetadata.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Return an object that calls 'describe' which can be streamed
 */
detail::DescribeMetadata
to_stream(const UniverseMetadata& md, const Tracker& tr)
{
    return detail::DescribeMetadata{md, tr};
}

//---------------------------------------------------------------------------//
/*!
 * Return a streamable object for metadata including daughter content.
 *
 * In addition to the basic output, this describes daughters and cells that are
 * local to the object.
 */
detail::DescribeUniverseMetadata
to_stream(const UniverseMetadata&                          md,
          const Universe&                                  u,
          std::function<const ObjectMetadata&(UniverseId)> get_daughter_md,
          std::function<void(std::ostream&, VolumeId)>     print_local_vol)
{
    CELER_EXPECT(get_daughter_md);
    CELER_EXPECT(print_local_vol);
    return detail::DescribeUniverseMetadata{
        md, u, std::move(get_daughter_md), std::move(print_local_vol)};
}

//---------------------------------------------------------------------------//
} // namespace celeritas
