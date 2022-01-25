//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/detail/DescribeMetadata.cc
 * \brief DescribeMetadata class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "DescribeMetadata.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include "base/Assert.hh"
#include "base/Range.hh"
#include "orange/track/Universe.hh"
#include "../ObjectMetadata.hh"
#include "../UniverseMetadata.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
std::ostream& operator<<(std::ostream& os, const DescribeMetadata& dm)
{
    dm.md.describe(os, dm.tracker);
    return os;
}

//---------------------------------------------------------------------------//
std::ostream& operator<<(std::ostream& os, const DescribeUniverseMetadata& dm)
{
    CELER_EXPECT(dm.get_daughter_md);
    CELER_EXPECT(dm.print_local_vol);
    CELER_EXPECT(dm.universe.num_volumes() > 0);

    // Write universe contents (mesh/surfaces/volumes)
    dm.md.describe(os, dm.universe.tracker());

    // Figure out cell lengths
    std::vector<std::string> id_to_labels(dm.universe.num_volumes());
    for (auto local_vol : range(dm.universe.num_volumes()))
    {
        id_to_labels[local_vol] = dm.md.id_to_label(VolumeId{local_vol});
    }
    size_type cell_width
        = std::max_element(id_to_labels.begin(),
                           id_to_labels.end(),
                           [](const std::string& lhs, const std::string& rhs) {
                               return lhs.size() < rhs.size();
                           })
              ->size();
    cell_width = std::min<size_type>(cell_width, 40);

    std::string daughter_sep(80, '=');
    daughter_sep[cell_width] = ' ';
    daughter_sep.back()      = '\n';

    // Write daughters
    os << daughter_sep << std::setw(cell_width) << "Cell"
       << " Fill\n"
       << daughter_sep;
    for (auto local_vol : range(dm.universe.num_volumes()))
    {
        os << std::setw(cell_width) << dm.md.id_to_label(VolumeId{local_vol})
           << ' ';
        const auto* d = dm.universe.daughter(VolumeId{local_vol});
        if (d)
        {
            // Print daughter name and transformation
            CELER_ASSERT(d->universe);
            os << dm.get_daughter_md(d->universe->id()).name();
            if (d->transform.is_nontrivial())
            {
                os << ' ' << d->transform;
            }
        }
        else
        {
            // Print cell contents
            dm.print_local_vol(os, VolumeId{local_vol});
        }
        os << '\n';
    }
    os << daughter_sep << '\n';

    return os;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
