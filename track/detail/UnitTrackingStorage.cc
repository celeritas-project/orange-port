//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/UnitTrackingStorage.cc
 * \brief UnitTrackingStorage class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "UnitTrackingStorage.hh"

#include <vector>
#include "Nemesis/comm/Timing.hh"
#include "base/Range.hh"
#include "orange/kdtree/KDBuilder.hh"
#include "orange/surfaces/SurfaceAction.hh"
#include "orange/Fuzziness.hh"
#include "orange/construct/UnitRegion.hh"
#include "SurfaceFunctors.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
UnitTrackingStorage::UnitTrackingStorage(SurfaceContainer surfaces,
                                         const std::vector<UnitRegion>& regions)
    : surfaces_(std::move(surfaces)), cells_(regions)
{
    CELER_EXPECT(cells_.num_surfaces() <= surfaces_.size());

    using geometria::KDBuilder;

    SCOPED_TIMER(
        "celeritas::detail::UnitTrackingStorage::UnitTrackingStorage");

    // Surface -> cell connectivity
    std::vector<std::vector<volume_int>> connectivity;
    connectivity.resize(cells_.num_surfaces());

    // Bounding box should be expanded to avoid connectivity problems due to
    // rounding point error
    KDBuilder::Params params;
    params.expand_bbox = fuzziness().bbox_rel();
    KDBuilder kdbuilder(std::move(params));

    // Functor to calculate number of surface intersections in a cell
    auto calc_num_intersections
        = make_surface_action(surfaces_, detail::NumIntersections());

    for (auto volume_idx : range(regions.size()))
    {
        const auto& rgn      = regions[volume_idx];
        const auto& cell_def = cells_[VolumeId(volume_idx)];

        // Add to kdtree if it's a region that can accept particles
        if (cell_def.zorder
            != static_cast<zorder_int>(ZOrder::implicit_exterior))
        {
            kdbuilder.insert(volume_idx, rgn.bbox);
        }

        // Update surface-related properties for cell
        size_type num_intersections = 0;
        for (SurfaceId sid : cell_def.faces)
        {
            // Increment number of intersections for this surface
            num_intersections += calc_num_intersections(sid);
            CELER_ASSERT(sid < connectivity.size());

            // Add connectivity from surface to the current cell
            connectivity[sid.get()].push_back(volume_idx);
        }
        CELER_ASSERT(num_intersections >= cell_def.faces.size());
        cells_.set_num_intersections(VolumeId(volume_idx), num_intersections);

        // Set flags
        CellContainer::flag_int flags = 0;
        if (rgn.has_internal_surfaces || cell_def.faces.empty())
        {
            // Internal surfaces if specified *or* for complex logic treatment
            // if it has no faces at all (infinite region)
            flags |= CellContainer::INTERNAL_SURFACES;
        }
        cells_.set_flags(VolumeId(volume_idx), flags);
    }

    // Flatten connectivity: since cell IDs were iterated in ascending order,
    // the per-surface connectivity is also in ascending order
    CELER_ASSERT(std::all_of(connectivity.begin(),
                             connectivity.end(),
                             [](const std::vector<volume_int>& vec_id) {
                                 return std::is_sorted(vec_id.begin(),
                                                       vec_id.end());
                             }));
    connectivity_ = TableCellId(connectivity.begin(), connectivity.end());

    // Complete kd-tree
    kdbuilder.build(kdtree_);
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
