//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/UnitMetadata.i.hh
 * \brief UnitMetadata inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Get detailed metadata for a surface
 */
auto UnitMetadata::surface_md(SurfaceId id) const -> SpanConstShapeFace
{
    CELER_EXPECT(id < md_.surfaces.size());
    auto view = md_.surfaces[id.get()];

    return {view.data(), view.size()};
}

//---------------------------------------------------------------------------//
/*!
 * Get detailed metadata for a cell
 */
const ObjectMetadata& UnitMetadata::vol_md(VolumeId id) const
{
    CELER_EXPECT(id < md_.cells.size());
    return md_.cells[id.get()];
}

//---------------------------------------------------------------------------//
} // namespace celeritas
