//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/GeometryParams.hh
 * \brief GeometryParams class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <map>
#include <memory>
#include <vector>
#include "orange/Definitions.hh"
#include "orange/query/UniverseMetadata.hh"
#include "orange/track/Universe.hh"
#include "../Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * \struct GeometryParams
 * Input parameters for constructing an ORANGE geometry.
 *
 * This is the final output of the ORANGE code before the Shift geometry
 * adapters are used. It is the output of the GeometryBuilder.
 */
//---------------------------------------------------------------------------//
struct GeometryParams
{
    using BoundaryState   = geometria::BoundaryState;
    using SPConstUniverse = std::shared_ptr<const Universe>;
    using SPMetadata      = std::shared_ptr<UniverseMetadata>;

    using MapSurfaceBoundary = std::map<SurfaceId, BoundaryState>;
    using VecUniverse        = std::vector<SPConstUniverse>;
    using VecMatid           = std::vector<geometria::matid_type>;
    using VecMetadata        = std::vector<SPMetadata>;

    VecUniverse        universes;
    VecMetadata        md;
    VecMatid           matids;
    MapSurfaceBoundary boundaries;

    size_type max_duplicate_cell_warn = 0;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
