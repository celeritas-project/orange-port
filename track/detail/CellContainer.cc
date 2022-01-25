//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/CellContainer.cc
 * \brief CellContainer class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CellContainer.hh"

#include "base/Range.hh"
#include "orange/construct/UnitRegion.hh"
#include "BuildLogic.hh"

using make_span;

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Construct with regions.
 *
 * This flattens the list of surface IDs per cell and converts the logic
 * array to use those cell indexes. Since integer sizes may be different
 * between the input surface ID type and the local "logic_int" storage we also
 * include checks for those limits.
 */
CellContainer::CellContainer(const std::vector<UnitRegion>& regions)
{
    CELER_EXPECT(!regions.empty());
    CELER_EXPECT(
        std::is_sorted(regions.begin(), regions.end(), ZorderPriority()));

    // Vector of # surfaces and # logic for each cell
    std::vector<std::pair<size_type, size_type>> sizes;

    for (const UnitRegion& rgn : regions)
    {
        CELER_ASSERT(rgn.interior);
        // Add surface IDs to list
        auto ids = rgn.interior.get_local_surfids();
        Insist(ids.size() < LogicEvaluator::max_num_surfaces(),
               "Too many surfaces (" << ids.size() << " >= "
                                     << LogicEvaluator::max_num_surfaces()
                                     << ") in a region");
        surfids_.insert(surfids_.end(), ids.begin(), ids.end());
        if (!ids.empty())
        {
            num_surfaces_ = std::max(num_surfaces_, ids.back().get() + 1);
        }

        // Convert user logic vector to internal logic
        const auto& rgn_logic       = rgn.interior.logic();
        auto        logic_start_idx = logic_.size();
        logic_.resize(logic_start_idx + rgn_logic.size());
        std::transform(rgn_logic.begin(),
                       rgn_logic.end(),
                       logic_.begin() + logic_start_idx,
                       BuildLogic{make_span(ids)});

        // Add insertions to vector of temporary sizes
        sizes.push_back({ids.size(), rgn_logic.size()});
    }
    CELER_ASSERT(sizes.size() == regions.size());

    // Allocate cells
    cells_.resize(regions.size());

    // Fill cell storage
    const SurfaceId* surf_start = surfids_.data();
    const logic_int* lgc_start  = logic_.data();
    for (auto i : range(regions.size()))
    {
        const auto& size_surf_lgc = sizes[i];
        const auto& rgn           = regions[i];
        auto&       local_vol     = cells_[i];

        // Only VolumeId{0} is allowed to be exterior
        CELER_ASSERT((rgn.zorder != ZOrder::exterior
                      && rgn.zorder != ZOrder::implicit_exterior)
                     || i == 0);

        // Create local cell
        local_vol.zorder = static_cast<zorder_int>(rgn.zorder);
        local_vol.faces  = {surf_start, size_surf_lgc.first};
        local_vol.logic  = {lgc_start, size_surf_lgc.second};

        // Update start pointers
        surf_start = local_vol.faces.end();
        lgc_start  = local_vol.logic.end();
    }

    CELER_ENSURE(surf_start == surfids_.data() + surfids_.size());
    CELER_ENSURE(lgc_start == logic_.data() + logic_.size());
}

//---------------------------------------------------------------------------//
/*!
 * Reconstruct a CSG cell for debug output
 */
UnitRegion CellContainer::get_region(VolumeId id) const
{
    CELER_EXPECT(id < cells_.size());

    // Input and output
    const mapped_type& cell = (*this)[id];
    UnitRegion         result;

    // Convert face IDs into surface IDs
    CSGCell::VecLogic logic(cell.logic.size());
    std::transform(cell.logic.begin(),
                   cell.logic.end(),
                   logic.begin(),
                   UnbuildLogic{cell.faces});
    result.interior = std::move(logic);

    // Copy flags
    result.zorder                = static_cast<ZOrder>(cell.zorder);
    result.has_internal_surfaces = cell.flags & INTERNAL_SURFACES;

    CELER_ENSURE(result.interior);
    return result;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
