//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/GeometryBuilder.hh
 * \brief GeometryBuilder class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <map>
#include <memory>
#include <unordered_map>
#include "base/FastHashMap.hh"
#include "orange/Definitions.hh"
#include "GeometryParams.hh"
#include "Proto.hh"

namespace celeritas
{
class PlacedShape;
class Universe;
//---------------------------------------------------------------------------//
/*!
 * Construct an ORANGE geometry class from protos.
 *
 * \code
    GeometryBuilder build_geometry;
    auto result = build_geometry(global_proto, boundary_conditions);
    return result.geometry;
 \endcode
 */
class GeometryBuilder
{
  public:
    //@{
    //! Public type aliases
    using SPConstProto         = std::shared_ptr<const Proto>;
    using SPConstShape         = std::shared_ptr<const PlacedShape>;
    using ShapeFace            = std::pair<SPConstShape, std::string>;
    using BoundaryState        = geometria::BoundaryState;
    using MapProtoUniv         = std::unordered_map<SPConstProto, UniverseId>;
    using MapShapeFaceBoundary = std::map<ShapeFace, BoundaryState>;
    //@}

    struct result_type
    {
        MapProtoUniv   proto_to_univ; //!< Proto/universe mapping
        GeometryParams geometry_params;
    };

  public:
    // Constructor takes options (no options yet? maybe tolerances?)
    GeometryBuilder() = default;

    // Build the geometry
    result_type
    operator()(SPConstProto global, MapShapeFaceBoundary boundaries) const;

  private:
    using VecBuiltProto      = std::vector<Proto::BuildResult>;
    using VecUniverse        = GeometryParams::VecUniverse;
    using VecMetadata        = GeometryParams::VecMetadata;
    using VecMatid           = GeometryParams::VecMatid;
    using MapSurfaceBoundary = GeometryParams::MapSurfaceBoundary;

    std::pair<MapProtoUniv, VecBuiltProto> build_protos(SPConstProto top) const;

    VecUniverse build_universes(const MapProtoUniv& proto_univ_ids,
                                VecBuiltProto&      built_protos) const;

    VecMetadata build_metadata(const MapProtoUniv& proto_univ_ids,
                               VecBuiltProto&      built_protos) const;

    VecMatid build_matids(const VecBuiltProto& built_protos,
                          const VecUniverse&   universes) const;

    MapSurfaceBoundary
    build_boundaries(const UniverseMetadata& top_md,
                     MapShapeFaceBoundary    face_boundaries) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
