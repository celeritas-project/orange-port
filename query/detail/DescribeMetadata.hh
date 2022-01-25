//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/detail/DescribeMetadata.hh
 * \brief DescribeMetadata class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <functional>
#include "orange/Definitions.hh"

namespace celeritas
{
class ObjectMetadata;
class UniverseMetadata;
class Tracker;
class Universe;
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Helper class for printing metadata and internals.
 */
struct DescribeMetadata
{
    const UniverseMetadata& md;
    const Tracker&          tracker;
};

//---------------------------------------------------------------------------//
/*!
 * Helper class that prints metadata and daughter information
 *
 * The \c print_not_daughter member will be called with the local cell ID of
 * any cell that does not have a daughter
 */
struct DescribeUniverseMetadata
{
    const UniverseMetadata&                          md;
    const Universe&                                  universe;
    std::function<const ObjectMetadata&(UniverseId)> get_daughter_md;
    std::function<void(std::ostream&, VolumeId)>     print_local_vol;
};

//---------------------------------------------------------------------------//

std::ostream& operator<<(std::ostream&, const DescribeMetadata&);
std::ostream& operator<<(std::ostream&, const DescribeUniverseMetadata&);

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
