//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/XMLBuilder.hh
 * \brief XMLBuilder class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "orange/ORANGEGeometry.hh"
#include "PlacedShape.hh"
#include "Proto.hh"
#include "UnitProto.hh"
#include "ArrayProto.hh"

namespace Teuchos
{
class ParameterList;
template<class T>
class Array;
} // namespace Teuchos

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Build an ORANGE geometry from an orangeinp XML file.
 */
class XMLBuilder
{
  public:
    //@{
    //! Public type aliases
    using ParameterList = Teuchos::ParameterList;
    using VecString     = std::vector<std::string>;
    using matid_type    = geometria::matid_type;
    using MapStrMatid   = std::unordered_map<std::string, matid_type>;
    using SPGeometry    = std::shared_ptr<ORANGEGeometry>;
    //@}

    //! Constructed result
    struct result_type
    {
        SPGeometry geo;   //!< Geometry
        VecString  comps; //!< Composition names
    };

  public:
    // Construct with names for the compositions, needed for Omnibus
    explicit XMLBuilder(MapStrMatid comp_names);

    // Construct with names where input names should be matids
    XMLBuilder();

    // Default destructor
    ~XMLBuilder();

    // Build from filename
    result_type operator()(const char* path);

    // Build from a parameter list
    result_type operator()(const ParameterList& plist);

  private:
    //// TYPES ////

    using ArrInt       = Teuchos::Array<int>;
    using ArrDbl       = Teuchos::Array<real_type>;
    using ArrStr       = Teuchos::Array<std::string>;
    using SPConstProto = std::shared_ptr<Proto>;
    using SPConstArray = std::shared_ptr<ArrayProto>;
    using SPConstShape = std::shared_ptr<const PlacedShape>;
    using ShapeFace    = std::pair<SPConstShape, std::string>;
    using VecHole      = std::vector<UnitProto::Hole>;
    using RegionVec    = Proto::RegionVec;
    using ProtoVec3    = ArrayProto::ProtoVec3;
    using MapStrMd     = std::unordered_map<std::string, ObjectMetadata>;
    using MapStrShape  = std::unordered_map<std::string, SPConstShape>;
    using MapStrProto  = std::unordered_map<std::string, SPConstProto>;

    //// PERSISTENT DATA ////

    MapStrMatid external_compnames_;

    //// STATE DATA (should be rolled into the operator()) ////

    // Map of externally derived compositions
    MapStrMatid compnames_;

    // Map of name -> proto-universe
    MapStrProto protos_;

    // List of reflecting boundaries
    std::vector<ShapeFace> reflecting_;

    // Map of missing comps -> metadata in which they first appeared
    MapStrMd missing_comps_;

#if 0
    // Reactor core proto constructed in VERA
    static VecProto& reactor_core_protos();
#endif
    //// IMPLEMENTATION METHODS ////

    void build_comp_names(const ParameterList& plist);
    void build_proto(const ParameterList& plist);

    // Proto-universe builders
    SPConstProto build_unit_proto(const ParameterList& plist);
    SPConstProto build_random_proto(const ParameterList& plist);
    SPConstProto build_rect_array_proto(const ParameterList& plist);
    SPConstProto build_hex_array_proto(const ParameterList& plist);
    SPConstProto build_dod_array_proto(const ParameterList& plist);
    SPConstProto build_rtk_proto(const ParameterList& plist);
    SPConstProto build_core_proto(const ParameterList& plist);

    // Get the material ID from a composition name
    matid_type
    comp_to_matid(const std::string& compname, const ObjectMetadata& context);
    // Build a region definition vector from a list of senses+shapes
    RegionVec build_region_vec(const MapStrShape&, const ArrStr&) const;
    // Construct local shapes
    MapStrShape build_shapes(const ParameterList& plist);
    // Build holes in the current universe
    VecHole build_holes(const ParameterList& plist,
                        ZOrder               zorder,
                        MapStrShape*         shapes) const;
    // Build units in an array, reordering from "natural" to ijk
    void
    fill_units(const ArrStr& fill, const char* desc, ProtoVec3* units) const;
    // Convert arrays to masked units if boundary or holes are given
    SPConstProto build_ggarray(SPConstArray, const ParameterList&);
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
