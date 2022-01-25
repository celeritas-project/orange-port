//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file input/XMLBuilder.cc
 * \brief XMLBuilder class definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "XMLBuilder.hh"

#include <set>
#include <unordered_set>

#include <TeuchosParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Nemesis/comm/Logger.hh"
#include "Nemesis/comm/TeuchosComm.hh"
#include "Nemesis/comm/Timing.hh"
#include "Nemesis/database/PlistUtilities.hh"
#include "Nemesis/io_utils/FileFunctions.hh"
#include "base/Join.hh"
#include "base/Make_RCP.hh"
#include "base/RegularIndexer.hh"

#include "orange/Fuzziness.hh"

#include "DodeArrayProto.hh"
#include "HexArrayProto.hh"
#include "IntersectionShape.hh"
#include "GeometryBuilder.hh"
#include "RandomProto.hh"
#include "RectArrayProto.hh"
#include "FaceNameCalculator.hh"

#include "detail/XMLUtils.hh"

using SetString = std::set<std::string>;
using DimVector = celeritas::ArrayProto::DimVector;

namespace celeritas
{
//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
namespace
{
//---------------------------------------------------------------------------//
//! Construct useful metadata object for a unit's boundary
ObjectMetadata build_boundary_md(const ObjectMetadata& unit_md)
{
    CELER_EXPECT(unit_md);
    ObjectMetadata::Params params;
    params.name       = "[EXTERIOR]";
    params.provenance = unit_md.provenance();
    return ObjectMetadata{std::move(params)};
}

template<class Maplike>
struct StreamableInvalid
{
    const SetString& actual;
    const char*      actual_desc;
    const Maplike&   expected;
    const char*      expected_desc;
};

template<class Maplike>
std::ostream& operator<<(std::ostream& os, const StreamableInvalid<Maplike>& sm)
{
    os << "Invalid " << sm.actual_desc << " in " << sm.expected_desc << ": "
       << join(sm.actual.begin(), sm.actual.end(), ", ")
       << "; valid choices are "
       << join(sm.expected.begin(),
               sm.expected.end(),
               ", ",
               [](const typename Maplike::value_type& kv) { return kv.first; });
    return os;
}

//---------------------------------------------------------------------------//
} // namespace

//---------------------------------------------------------------------------//
/*!
 * Construct with predefined composition names.
 */
XMLBuilder::XMLBuilder(MapStrMatid comp_names)
    : external_compnames_(std::move(comp_names))
{
    std::vector<MapStrMatid::value_type> duplicates;
    std::unordered_set<matid_type>       ids;
    for (const auto& kv : external_compnames_)
    {
        auto iter_inserted = ids.insert(kv.second);
        if (!iter_inserted.second)
        {
            duplicates.push_back(kv);
        }
    }
    Validate(duplicates.empty(),
             "ORANGE XML compositions cannot share a matid: duplicates are "
                 << join_stream(
                        duplicates.begin(),
                        duplicates.end(),
                        ", ",
                        [](std::ostream& os, const MapStrMatid::value_type& kv) {
                            os << kv.second << " -> '" << kv.first << '\'';
                        }));
}

//---------------------------------------------------------------------------//
/*!
 * Default constructor defines no compositions
 */
XMLBuilder::XMLBuilder() = default;

//---------------------------------------------------------------------------//
//! Default destructor
XMLBuilder::~XMLBuilder() = default;

//---------------------------------------------------------------------------//
/*!
 * Build the geometry from an XML input file.
 */
auto XMLBuilder::operator()(const char* filename) -> result_type
{
    Teuchos::RCP<ParameterList> plist = make_rcp<ParameterList>("geometry");

    // Comm pointer (use the local MPI comm if possible)
    auto comm = TeuchosComm::get_default();
    CELER_ASSERT(!comm.is_null());

    if (comm->getRank() == 0)
    {
        CELER_VALIDATE(is_file(filename),
                       << "No ORANGE XML file is present at \"" << filename
                       << "\"");
    }

    // Load parameters from disk on processor zero and broadcast them
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        filename, plist.ptr(), *comm);
    CELER_ASSERT(!plist.is_null());

    return (*this)(*plist);
}

//---------------------------------------------------------------------------//
/*!
 * Build the geometry from a parameter list.
 */
auto XMLBuilder::operator()(const ParameterList& plist) -> result_type
{
    missing_comps_.clear();

    SCOPED_TIMER("celeritas::XMLBuilder::operator()");

    // Load global settings from the GEOMETRY list
    const ParameterList& geo_plist = plist.sublist("GEOMETRY");

    // Set options
    if (geo_plist.get<bool>("write_kdtree"))
    {
        log(WARNING) << "KD tree output is currently "
                        "unimplemented";
    }
    fuzziness() = Fuzziness(geo_plist.get<real_type>("tolerance"),
                            geo_plist.get<real_type>("length_scale"));

    // Load user-provided compositions if they exist
    compnames_ = external_compnames_;
    if (geo_plist.isParameter("composition"))
    {
        this->build_comp_names(geo_plist);
    }

    // Build universes
    for (const auto& kv : plist.sublist("UNIVERSE"))
    {
        this->build_proto(Teuchos::getValue<ParameterList>(kv.second));
    }

    // Validate that no comps were missing
    CELER_VALIDATE(
        missing_comps_.empty(),
        << "The following compositions do not exist (but were requested in "
           "at least one cell): "
        << join_stream(missing_comps_.begin(),
                       missing_comps_.end(),
                       ", ",
                       [](std::ostream& os, const MapStrMd::value_type& kv) {
                           os << kv.first << " (from " << kv.second << ")";
                       }));

    // Build geometry
    result_type result;
    {
        GeometryBuilder build_geo;

        // Find global universe
        const std::string& name = geo_plist.get<std::string>("global");
        auto               global_proto_iter = protos_.find(name);
        CELER_VALIDATE(global_proto_iter != protos_.end(),
                       << "Global universe '" << name << "' does not exist");

        // Construct boundaries
        GeometryBuilder::MapShapeFaceBoundary reflecting;
        for (auto& shape_face : reflecting_)
        {
            reflecting.emplace(std::move(shape_face), geometria::REFLECT);
        }

        // Build geometry
        auto built
            = build_geo(global_proto_iter->second, std::move(reflecting));
        result.geo = std::make_shared<ORANGEGeometry>(
            std::move(built.geometry_params));
    }

    // Build composition names
    result.comps.reserve(compnames_.size());
    for (const auto& kv : compnames_)
    {
        size_type matid = kv.second;
        if (matid >= result.comps.size())
        {
            result.comps.resize(matid + 1);
        }
        result.comps[matid] = kv.first;
    }

    // Clean up temporary state
    protos_.clear();
    reflecting_.clear();

    // Construct final universe
    CELER_ENSURE(result.geo);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Build compositions.
 */
void XMLBuilder::build_comp_names(const ParameterList& plist)
{
    CELER_EXPECT(plist.isParameter("composition"));
    CELER_EXPECT(plist.isParameter("matid"));

    const ArrStr& comp_names = plist.get<ArrStr>("composition");
    const ArrInt& matids     = plist.get<ArrInt>("matid");

    struct Dupe
    {
        std::string name;
        matid_type  existing_matid;
        matid_type  given_matid;
    };
    std::vector<Dupe> duplicates;

    auto matid_iter = matids.begin();
    for (const std::string& name : comp_names)
    {
        matid_type matid = *matid_iter++;
        auto       added = compnames_.insert({name, matid});
        if (!added.second && added.first->second != matid)
        {
            // Conflicting comp name
            duplicates.push_back(
                {added.first->first, added.first->second, matid});
        }
    }

    if (!duplicates.empty())
    {
        log(WARNING)
            << "Duplicate composition names were given in ORANGE input: "
            << join_stream(duplicates.begin(),
                           duplicates.end(),
                           ", ",
                           [](std::ostream& os, const Dupe& d) {
                               os << d.name << " -> {" << d.existing_matid
                                  << ',' << d.given_matid << '}';
                           });
    }
}

//---------------------------------------------------------------------------//
/*!
 * Build a single proto from the parameter list.
 */
void XMLBuilder::build_proto(const ParameterList& plist)
{
    using ProtoBuilder = SPConstProto (XMLBuilder::*)(const ParameterList&);
    using MapBuilders  = std::unordered_map<std::string, ProtoBuilder>;

    static const MapBuilders proto_builders = {
        // clang-format off
        {"unit",        &XMLBuilder::build_unit_proto},
        {"random",      &XMLBuilder::build_random_proto},
        {"array",       &XMLBuilder::build_rect_array_proto},
        {"hexarray",    &XMLBuilder::build_hex_array_proto},
        {"dodarray",    &XMLBuilder::build_dod_array_proto},
        {"rtk",         &XMLBuilder::build_rtk_proto},
        {"core",        &XMLBuilder::build_core_proto},
        // clang-format on
    };

    // Build the universe: use the "_type" parameter to determine what
    // function to call.
    auto iter = proto_builders.find(plist.get<std::string>("_type"));
    CELER_VALIDATE(iter != proto_builders.end(),
                   << "Invalid universe type '"
                   << plist.get<std::string>("_type") << "'");

    // Call member function pointer for the proto building method
    ProtoBuilder fp = iter->second;
    CELER_ASSERT(iter->second);
    SPConstProto proto = (this->*fp)(plist);
    CELER_ASSERT(proto);

    // Add proto-universe to our list
    auto result = protos_.insert({proto->metadata().name(), proto});
    CELER_VALIDATE(result.second,
                   << "Universe '" << proto->metadata().name()
                   << "' cannot have the same name as another "
                      "universe");
}

//---------------------------------------------------------------------------//
/*!
 * Build a "unit" proto.
 */
auto XMLBuilder::build_unit_proto(const ParameterList& plist) -> SPConstProto
{
    CELER_EXPECT(plist.get<std::string>("_type") == "unit");

    // Create universe
    UnitProto::Params params;
    params.md                 = detail::build_md(plist);
    params.otf_error_checking = plist.get<bool>("otf_error_checking");

    // Build shapes
    MapStrShape shapes = this->build_shapes(plist);

    // Build holes, save their shapes, update their zorder
    params.holes = this->build_holes(
        plist,
        plist.get<bool>("implicit_holes") ? ZOrder::hole : ZOrder::media,
        &shapes);

    // TODO: support placed arrays

    // Add cells (must be after both shapes and holes, because cells can refer
    // to both). It's optional because the user might want this universe to
    // consist of an entire other universe (e.g. if you want an 'array' to be
    // the external universe)
    if (plist.isSublist("CELL"))
    {
        for (const auto& kv : plist.sublist("CELL"))
        {
            const auto& plist = Teuchos::getValue<ParameterList>(kv.second);

            // Construct single "media" entry with metadata
            UnitProto::Media m;
            m.md = detail::build_md(plist);

            // Load surfaces
            m.interior
                = this->build_region_vec(shapes, plist.get<ArrStr>("shapes"));
            CELER_VALIDATE(!m.interior.empty(),
                           << "No surfaces were defined for cell " << m.md);

            // Transform composition to matid
            m.matid = this->comp_to_matid(
                plist.get<std::string>("composition"), m.md);

            // Get cell volume if provided
            if (plist.isParameter("volume"))
            {
                m.volume = plist.get<real_type>("volume");
            }

            params.media.push_back(std::move(m));
        }
    }

    // Build boundary definition
    CELER_VALIDATE(plist.isParameter("interior"),
                   << "Unit universe " << params.md
                   << " must have an 'interior' entry.");
    params.boundary.interior
        = this->build_region_vec(shapes, plist.get<ArrStr>("interior"));
    params.boundary.implicit_boundary = plist.get<bool>("implicit_boundary");
    params.boundary.md                = build_boundary_md(params.md);

    return std::make_shared<UnitProto>(std::move(params));
}

//---------------------------------------------------------------------------//
/*!
 * Build a 'random' universe of particles.
 */
auto XMLBuilder::build_random_proto(const ParameterList& plist) -> SPConstProto
{
    CELER_EXPECT(plist.get<std::string>("_type") == "random");

    // Create universe
    RandomProto::Params params;
    params.md = detail::build_md(plist);

    // Build options
    {
        const ParameterList& opt_list = plist.sublist("OPTIONS");
        RandomProto::Options options;
        options.seed               = opt_list.get<int>("seed");
        options.failure_batch_size = opt_list.get<int>("failure_batch_size");
        options.failure_tolerance
            = opt_list.get<real_type>("failure_tolerance");

        params.options = std::move(options);
    }

    // Build fill material
    params.fill_matid = this->comp_to_matid(
        plist.get<std::string>("composition"), params.md);

    // Build boundary shape
    const auto& shape_plist = plist.get<ParameterList>("SHAPE");
    params.interior = {{Sense::inside, detail::build_shape(shape_plist)}};

    // Add particle universes
    const ArrStr& universes = plist.get<ArrStr>("fill");
    const ArrDbl& vol_fracs = plist.get<ArrDbl>("volume_fraction");
    Insist(universes.size() == vol_fracs.size(),
           "Size mismatch between 'fill' and 'volume_fraction'");

    SetString missing;
    for (auto i : range(universes.size()))
    {
        auto iter = protos_.find(universes[i]);
        if (iter == protos_.end())
        {
            missing.insert(universes[i]);
            continue;
        }

        RandomProto::Particle p;
        p.proto           = iter->second;
        p.volume_fraction = vol_fracs[i];
        p.md              = iter->second->metadata();
        params.particles.push_back(std::move(p));
    }

    CELER_VALIDATE(missing.empty(),
                   << (StreamableInvalid<MapStrProto>{
                          missing, "daughters", protos_, "random universe"}));

    return std::make_shared<RandomProto>(std::move(params));
}

//---------------------------------------------------------------------------//
/*!
 * Build an array universe.
 */
auto XMLBuilder::build_rect_array_proto(const ParameterList& plist)
    -> SPConstProto
{
    CELER_EXPECT(plist.get<std::string>("_type") == "array");

    RectArrayProto::Params params;
    params.md = detail::build_md(plist);
    {
        // For embedding in GG universe
        ObjectMetadata::Params md;
        md.name        = params.md.name() + "+";
        md.description = "Embedded unplaced array";
        md.provenance  = params.md.provenance();
        params.md      = ObjectMetadata{std::move(md)};
    }

    // Build units
    params.units.resize(DimVector(
        plist.get<int>("nx"), plist.get<int>("ny"), plist.get<int>("nz")));
    this->fill_units(plist.get<ArrStr>("fill"), "rect array", &params.units);

    // Place in embedding Unit Proto
    return this->build_ggarray(
        std::make_shared<RectArrayProto>(std::move(params)), plist);
}

//---------------------------------------------------------------------------//
/*!
 * Build a hexagonal array universe.
 */
auto XMLBuilder::build_hex_array_proto(const ParameterList& plist)
    -> SPConstProto
{
    CELER_EXPECT(plist.get<std::string>("_type") == "hexarray");

    HexArrayProto::Params params;
    params.md = detail::build_md(plist);
    {
        // For embedding in GG universe
        ObjectMetadata::Params md;
        md.name        = params.md.name() + "+";
        md.description = "Embedded unplaced array";
        md.provenance  = params.md.provenance();
        params.md      = ObjectMetadata{std::move(md)};
    }

    const auto& layout_str = plist.get<std::string>("layout");
    if (layout_str == "rhomb")
    {
        params.layout = HexArrayProto::Layout::rhomboidal;
    }
    else if (layout_str == "rect")
    {
        params.layout = HexArrayProto::Layout::rectangular;
    }
    else
    {
        CELER_VALIDATE(false, << "Invalid hex layout " << layout_str);
    }

    // Build units
    params.units.resize(DimVector(
        plist.get<int>("nu"), plist.get<int>("nv"), plist.get<int>("nz")));
    this->fill_units(plist.get<ArrStr>("fill"), "hex array", &params.units);

    // Place in embedding Unit Proto
    return this->build_ggarray(
        std::make_shared<HexArrayProto>(std::move(params)), plist);
}

//---------------------------------------------------------------------------//
/*!
 * Build a rhombic dodecahedron array.
 */
auto XMLBuilder::build_dod_array_proto(const ParameterList& plist)
    -> SPConstProto
{
    CELER_EXPECT(plist.get<std::string>("_type") == "dodarray");

    DodeArrayProto::Params params;
    params.md = detail::build_md(plist);
    {
        // For embedding in GG universe
        ObjectMetadata::Params md;
        md.name        = params.md.name() + "+";
        md.description = "Embedded unplaced array";
        md.provenance  = params.md.provenance();
        params.md      = ObjectMetadata{std::move(md)};
    }

    // Build units
    params.units.resize(DimVector(
        plist.get<int>("nx"), plist.get<int>("ny"), plist.get<int>("nz")));
    this->fill_units(plist.get<ArrStr>("fill"), "dode array", &params.units);

    // Place in embedding Unit Proto
    return this->build_ggarray(
        std::make_shared<DodeArrayProto>(std::move(params)), plist);
}

//---------------------------------------------------------------------------//
/*!
 * Build an RTK array universe.
 */
auto XMLBuilder::build_rtk_proto(const ParameterList& plist) -> SPConstProto
{
    CELER_EXPECT(plist.get<std::string>("_type") == "rtk");

    NotImplemented("RTK proto for ORANGE");
}

//---------------------------------------------------------------------------//
/*!
 * Insert an VERA-built reactor for ex core dosimetry.
 */
auto XMLBuilder::build_core_proto(const ParameterList& plist) -> SPConstProto
{
    NotImplemented("VERA core construction for ORANGE");
}

//---------------------------------------------------------------------------//
/*!
 * Get the material ID from a composition name.
 */
auto XMLBuilder::comp_to_matid(const std::string&    compname,
                               const ObjectMetadata& context) -> matid_type
{
    matid_type matid = geometria::invalid_matid();

    if (!compnames_.empty())
    {
        auto iter = compnames_.find(compname);
        if (iter != compnames_.end())
        {
            matid = iter->second;
        }
        else
        {
            // Note that the comp is missing, validate at end of build
            missing_comps_.insert({compname, context});
        }
    }
    else
    {
        try
        {
            matid = std::stoi(compname);
        }
        catch (const std::invalid_argument& e)
        {
            missing_comps_.insert({compname, context});
            log(WARNING)
                << "Could not convert 'composition' entry '" << compname
                << "' in " << context
                << " to a numerical matid (and no named compositions exist)";
        }
    }
    return matid;
}

//---------------------------------------------------------------------------//
/*!
 * Build a region definition vector from a list of senses+shapes.
 *
 * The input processor compresses shapes plus senses into shapes with a
 * trailing `~`.
 */
auto XMLBuilder::build_region_vec(const MapStrShape& shapes,
                                  const ArrStr&      rdv) const -> RegionVec
{
    RegionVec result;
    result.reserve(rdv.size());

    SetString missing;
    for (std::string name : rdv)
    {
        CELER_VALIDATE(!name.empty(), << "Empty name for shape in region vec");

        // Extract sense
        Sense sense = Sense::inside;
        if (name.back() == '~')
        {
            sense = Sense::outside;
            name.pop_back();
        }

        // Extract shape ID
        auto iter = shapes.find(name);
        if (iter == shapes.end())
        {
            missing.insert(std::move(name));
            continue;
        }

        // Add to surface list
        result.push_back({sense, iter->second});
    }

    CELER_VALIDATE(
        missing.empty(),
        << (StreamableInvalid<MapStrShape>{
               missing, "shapes", shapes, "region definition vector"}));
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Construct shapes from the parent plist, or return empty if none.
 */
auto XMLBuilder::build_shapes(const ParameterList& plist) -> MapStrShape
{
    // Add shapes
    if (!plist.isSublist("SHAPE"))
    {
        return {};
    }

    MapStrShape shapes;
    for (const auto& kv : plist.sublist("SHAPE"))
    {
        const ParameterList& shape_plist
            = Teuchos::getValue<ParameterList>(kv.second);
        auto shape = detail::build_shape(shape_plist);
        CELER_ASSERT(shape);

        // Add any user-specified boundaries for the shape
        auto refl = get_default(shape_plist, "reflect", VecString{});
        if (refl.size() == 1 && refl.front() == "*")
        {
            // Add all faces
            refl = FaceNameCalculator()(*shape->shape());
            CELER_ASSERT(!refl.empty());
        }
        for (const std::string& surf : refl)
        {
            reflecting_.push_back({shape, surf});
        }

        // Store the shape ID added
        auto result = shapes.insert({shape->name(), std::move(shape)});
        CELER_VALIDATE(result.second,
                       << "Shape " << shape->metadata()
                       << " cannot have the same name as another "
                          "shape in the current universe");
    }
    return shapes;
}

//---------------------------------------------------------------------------//
/*!
 * Build holes that are referenceable as shapes.
 */
auto XMLBuilder::build_holes(const ParameterList& plist,
                             ZOrder               zorder,
                             MapStrShape*         shapes) const -> VecHole
{
    if (!plist.isSublist("HOLE"))
    {
        return {};
    }

    VecHole result;

    for (const auto& kv : plist.sublist("HOLE"))
    {
        const auto& plist = Teuchos::getValue<ParameterList>(kv.second);

        UnitProto::Hole h;
        h.md = detail::build_md(plist);

        // Extract fill proto
        const std::string& fill_name = plist.get<std::string>("fill");
        auto               iter      = protos_.find(fill_name);
        CELER_VALIDATE(iter != protos_.end(),
                       << "Given fill '" << fill_name
                       << "' must be an existing universe");
        h.proto = iter->second;

        // Apply transform and set zorder
        h.transform = detail::build_transform(plist);
        h.zorder    = zorder;

        // Add shape
        if (shapes)
        {
            auto result
                = shapes->insert({h.md.name(), UnitProto::make_hole_shape(h)});
            CELER_VALIDATE(result.second,
                           << "Hole " << h.md
                           << " cannot have the same name as another shape "
                              "in the current universe");
        }

        CELER_ASSERT(h);
        result.push_back(std::move(h));
    }

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Build units in an array, reordering from "natural" to ijk.
 *
 * This loops over elements: we want the input to be laid out in text like an
 * image, earlier in the document is higher K planes (slowest moving), then
 * earlier in a plane is a higher J row, then along a line is higher I row
 * (fastest moving). So we need to flip K and J indices from top to bottom, and
 * reorder the input KJI -> to destination IJK.
 *
 * This also checks for missing entries.
 */
void XMLBuilder::fill_units(const ArrStr& fill,
                            const char*   arr_type,
                            ProtoVec3*    units) const
{
    CELER_EXPECT(units);

    using def::I;
    using def::J;
    using def::K;

    // Index into the input list
    const auto&                              dims = units->dims();
    RegularIndexer<DimVector::value_type, 3> list_indexer(
        {dims[K], dims[J], dims[I]});
    DimVector list_index;
    CELER_ASSERT(list_indexer.size() == fill.size());

    SetString missing;

    for (auto k : range(dims[K]))
    {
        list_index[0] = list_indexer.dims()[0] - k - 1;
        for (auto j : range(dims[J]))
        {
            list_index[1] = list_indexer.dims()[1] - j - 1;
            for (auto i : range(dims[I]))
            {
                list_index[2] = i;

                // Extract fill name at the location in the list
                auto               list_idx  = list_indexer.index(list_index);
                const std::string& fill_name = fill[list_idx];

                // Get proto-universe
                auto iter = protos_.find(fill_name);
                if (iter == protos_.end())
                {
                    missing.insert(fill_name);
                    continue;
                }

                // Fill the array cell
                (*units)[{i, j, k}] = iter->second;
            }
        }
    }

    CELER_VALIDATE(missing.empty(),
                   << (StreamableInvalid<MapStrProto>{
                          missing, "units", protos_, arr_type}));
}

//---------------------------------------------------------------------------//
/*!
 * Transform an array into a unit enclosing an array.
 *
 * This is for backward compatibility with GG. We should allow "placed" arrays
 * later for optimization.
 */
auto XMLBuilder::build_ggarray(SPConstArray         arr_proto,
                               const ParameterList& plist) -> SPConstProto
{
    CELER_EXPECT(arr_proto);

    UnitProto::Params params;
    params.md = detail::build_md(plist);

    // Construct potential hole/boundary shapes
    MapStrShape shapes = this->build_shapes(plist);

    // Build holes (implicitly masking the array contents)
    params.holes = this->build_holes(plist, ZOrder::hole, &shapes);

    // Construct array
    {
        UnitProto::Array arr;
        arr.proto = arr_proto;

        // If array placement is 'unit', translate the array so that the
        // specified origin is the origin of the unit
        const std::string& origin_spec = plist.get<std::string>("origin_is");
        const Real3        origin = detail::get_space_vector(plist, "origin");
        if (origin_spec == "array")
        {
            // Lower left corner
            arr.transform = arr_proto->calc_placement(origin);
        }
        else if (origin_spec == "unit")
        {
            arr.transform = arr_proto->calc_placement(
                detail::get_dim_vector(plist, "place"), origin);
        }
        else
        {
            CELER_VALIDATE(
                false, << "Invalid value for 'origin_is': " << origin_spec);
        }

        if (plist.isParameter("interior"))
        {
            RegionVec interior = this->build_region_vec(
                shapes, plist.get<ArrStr>("interior"));
            ObjectMetadata interior_md = params.md;
            if (interior.size() == 1)
            {
                // Use metadata from single bounding shape
                interior_md = interior.front().second->metadata();
            }
            // Don't transform explicit shape
            auto shape = shape_from_rdv(std::move(interior), {}, interior_md);
            arr.interior = {{neg, shape}};
        }
        else
        {
            // Implicitly create boundary from array daughter
            auto shape = shape_from_rdv(
                arr_proto->interior(), arr.transform, params.md);
            arr.interior = {{neg, shape}};
        }
        arr.md = arr_proto->metadata();
        params.arrays.push_back(std::move(arr));
    }

    // Build boundary from array boundary
    params.boundary.interior          = params.arrays.back().interior;
    params.boundary.implicit_boundary = false;
    params.boundary.md                = build_boundary_md(params.md);

    return std::make_shared<UnitProto>(std::move(params));
}

//---------------------------------------------------------------------------//
} // namespace celeritas
