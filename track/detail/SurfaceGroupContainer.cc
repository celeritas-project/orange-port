//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/SurfaceGroupContainer.cc
 * \brief SurfaceGroupContainer class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SurfaceGroupContainer.hh"

#include <map>
#include "base/FastHashSet.hh"
#include "base/Range.hh"
#include "orange/surfaces/SurfaceAction.hh"
#include "SurfaceFunctors.hh"
#include "UnitTrackingStorage.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Construct with regions
 */
SurfaceGroupContainer::SurfaceGroupContainer(
    const UnitTrackingStorage& unit_storage)
{
    // Build map of
    using SetSurface       = FastHashSet<SurfaceId>;
    using DecreasingZorder = std::greater<zorder_int>;
    using MapZorderSurfs = std::map<zorder_int, SetSurface, DecreasingZorder>;
    MapZorderSurfs surfaces_by_zorder;
    for (auto c : range(unit_storage.cells().size()))
    {
        const auto& cell_def = unit_storage.cells()[VolumeId{c}];
        auto&       surfaces = surfaces_by_zorder[cell_def.zorder];
        for (SurfaceId sid : cell_def.faces)
        {
            surfaces.insert(sid);
        }
    }

    // Resize/reserve data
    groups_.reserve(unit_storage.surfaces().size());
    surface_ids_.resize(std::accumulate(
        surfaces_by_zorder.begin(),
        surfaces_by_zorder.end(),
        size_type(0),
        [](size_type size, const MapZorderSurfs::value_type& kv) {
            return size + kv.second.size();
        }));

    SurfaceId* surf_start             = surface_ids_.data();
    auto       calc_num_intersections = make_surface_action(
        unit_storage.surfaces(), detail::NumIntersections());

    // Add surfaces in decreasing zorder
    for (const auto& kv : surfaces_by_zorder)
    {
        // Temporary mutable view to the surfaces
        span<SurfaceId> faces{surf_start, kv.second.size()};
        surf_start += faces.size();

        // Append and sort the surfaces for this zorder
        std::copy(kv.second.begin(), kv.second.end(), faces.begin());
        std::sort(faces.begin(), faces.end());

        // Fill metadata and calculate number of intersections for this zorder
        value_type surfaces;
        surfaces.zorder            = kv.first;
        surfaces.faces             = faces;
        surfaces.num_intersections = 0;
        for (SurfaceId sid : surfaces.faces)
        {
            // Increment number of intersections for this surface
            surfaces.num_intersections += calc_num_intersections(sid);
        }

        groups_.push_back(surfaces);
    }
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
