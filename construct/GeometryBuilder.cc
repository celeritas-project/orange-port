//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/GeometryBuilder.cc
 * \brief GeometryBuilder class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "GeometryBuilder.hh"

#include <deque>
#include <iostream>
#include "base/FastHashSet.hh"
#include "base/Join.hh"
#include "base/Range.hh"
#include "orange/query/UnitMetadata.hh"
#include "orange/track/Universe.hh"
#include "PlacedShape.hh"
#include "Proto.hh"

namespace
{
//---------------------------------------------------------------------------//
using ShapeFace = celeritas::GeometryBuilder::ShapeFace;

//---------------------------------------------------------------------------//
/*!
 * Cast an item to bool, for use in \c std::all_of
 */
struct is_true
{
    template<class T>
    bool operator()(const T& obj)
    {
        return static_cast<bool>(obj);
    }
};

//---------------------------------------------------------------------------//
/*!
 * \brief Print a shape/face pair to stream.
 */
std::string to_string(const ShapeFace& face)
{
    CELER_EXPECT(face.first);
    std::string result(face.first->name());
    result += '.';
    result += face.second;
    return result;
}

//---------------------------------------------------------------------------//
} // namespace

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Build a geometry
 *
 * - Breadth-first construction of protos into \c Proto::BuildResult, which
 *   contains vector of daughter protos
 * - Universe IDs correspond to breadth-first order
 * - Reversed depth-first construction of universes (so that all daughters are
 *   guaranteed to have been created)
 */
auto GeometryBuilder::operator()(SPConstProto         top,
                                 MapShapeFaceBoundary boundaries) const
    -> result_type
{
    CELER_EXPECT(top);

    MapProtoUniv proto_univ_ids;

    // Build trackers and metadata from the protos
    VecBuiltProto built_protos;
    std::tie(proto_univ_ids, built_protos) = this->build_protos(std::move(top));
    CELER_ASSERT(proto_univ_ids.size() == built_protos.size());

    // Connect the layers via universes
    GeometryParams geo_params;
    geo_params.universes = this->build_universes(proto_univ_ids, built_protos);
    CELER_ASSERT(geo_params.universes.size() == proto_univ_ids.size());

    // Build metadata
    geo_params.md = this->build_metadata(proto_univ_ids, built_protos);
    CELER_ASSERT(geo_params.md.size() == proto_univ_ids.size());

    // Build materials
    geo_params.matids = this->build_matids(built_protos, geo_params.universes);

    // Build boundaries
    CELER_ASSERT(!geo_params.universes.empty() && geo_params.md.front());
    geo_params.boundaries = this->build_boundaries(*geo_params.md.front(),
                                                   std::move(boundaries));

    //// CONSTRUCT GEOMETRY ////

    return {std::move(proto_univ_ids), std::move(geo_params)};
}

//---------------------------------------------------------------------------//
/*!
 * Build protos, constructing trackers, metadata, and universe IDs.
 */
auto GeometryBuilder::build_protos(SPConstProto top) const
    -> std::pair<MapProtoUniv, VecBuiltProto>
{
    // Queue of unbuilt protos
    // TODO: make this a queue of vectors so that we can print out better
    // failure messages while building
    std::deque<const Proto*> queue{top.get()};
    // Map of universe IDs for each proto (top is zero)
    MapProtoUniv proto_univ_ids{{std::move(top), UniverseId{0}}};
    // Result of proto "build" corresponding to each universe ID
    VecBuiltProto build_results;

    // Build arguments for the next proto
    Proto::BuildArgs args;
    args.implicit_boundary = false;

    while (!queue.empty())
    {
        build_results.push_back(queue.front()->build(args));
        queue.pop_front();

        for (const auto& cell_daughter : build_results.back().daughters)
        {
            const SPConstProto& daughter_proto = cell_daughter.second.proto;
            CELER_ASSERT(daughter_proto);
            // Try inserting the next universe ID for this proto
            auto iter_inserted = proto_univ_ids.insert(
                {daughter_proto, UniverseId(proto_univ_ids.size())});
            if (iter_inserted.second)
            {
                // Insertion success: the daughter hasn't been visited, so add
                // it.
                queue.push_back(daughter_proto.get());
            }
        }

        // Allow higher universe to truncate all but top universe
        args.implicit_boundary = true;
    }
    CELER_ENSURE(build_results.size() == proto_univ_ids.size());

    return {std::move(proto_univ_ids), std::move(build_results)};
}

//---------------------------------------------------------------------------//
/*!
 * Construct universes from built trackers.
 *
 * The "built protos" use breadth-first traversal but to construct universes we
 * need to do reverse depth-first (to guarantee that daughters are built before
 * their parents, which need to reference them). Construct an ordering that
 * guarantees daughters will be built before parents, then build based on that
 * order.
 */
auto GeometryBuilder::build_universes(const MapProtoUniv& proto_univ_ids,
                                      VecBuiltProto&      build_results) const
    -> VecUniverse
{
    CELER_EXPECT(!proto_univ_ids.empty());
    CELER_EXPECT(build_results.size() == proto_univ_ids.size());

    using size_type = UniverseId::size_type;

    // Ordering is the *reverse* of a topological ordering of the proto
    // dependency DAG.
    std::vector<size_type> ordering;
    {
        // Set of protos visited and indices added
        FastHashSet<const Proto*> visited;
        FastHashSet<size_type>    added;
        visited.reserve(build_results.size());
        added.reserve(build_results.size());

        // Stack of daughters being visited: start with head node zero
        std::deque<size_type> stack{0};
        while (!stack.empty())
        {
            auto parent_id = stack.back();
            if (added.count(parent_id))
            {
                // Append to the ordering (reverse depth-first)
                ordering.push_back(parent_id);
                stack.pop_back();
                continue;
            }

            bool all_daughters_built = true;

            // Loop over daughters from proto build result
            for (const auto& cell_daughter : build_results[parent_id].daughters)
            {
                // Find universe ID corresponding to daughter
                const auto& sp_daughter    = cell_daughter.second.proto;
                auto        iter_insertion = visited.insert(sp_daughter.get());
                if (iter_insertion.second)
                {
                    // Add daughter ID to stack
                    auto daughter_id = proto_univ_ids.find(sp_daughter);
                    CELER_ASSERT(daughter_id != proto_univ_ids.end());
                    stack.push_back(daughter_id->second.get());
                    all_daughters_built = false;
                }
            }

            if (all_daughters_built && !added.count(parent_id))
            {
                added.insert(parent_id);
            }
        }

        CELER_ASSERT(added.size() == build_results.size());
        CELER_ASSERT(ordering.size() == build_results.size());
        CELER_ASSERT(
            FastHashSet<size_type>(ordering.begin(), ordering.end()).size()
            == ordering.size());
    }

    VecUniverse universes(build_results.size());

    // Loop over all universes in the constructed order
    for (size_type parent_idx : ordering)
    {
        auto& built = build_results[parent_idx];
        CELER_ASSERT(built.tracker);

        // Convert proto build result into universe build arguments by mapping
        // daughter protos to universe IDs to universe
        Universe::Params params;
        params.id      = UniverseId(parent_idx);
        params.tracker = std::move(built.tracker);

        if (!built.daughters.empty())
        {
            params.daughters.resize(params.tracker->num_volumes());
        }

        // Loop over cell/proto/transform from proto build result
        for (const auto& cell_daughter : built.daughters)
        {
            CELER_ASSERT(cell_daughter.first < params.daughters.size());
            // Find the universe ID corresponding to the daughter proto
            auto daughter_id = proto_univ_ids.find(cell_daughter.second.proto);
            CELER_ASSERT(daughter_id != proto_univ_ids.end());
            auto daughter_idx = daughter_id->second.get();

            // Insert daughter universe into the new parameters
            Universe::Daughter new_daughter;
            new_daughter.universe = universes[daughter_idx];
            CELER_ASSERT(new_daughter.universe);
            new_daughter.transform = cell_daughter.second.transform;
            params.daughters[cell_daughter.first.unchecked_get()]
                = std::move(new_daughter);
        }

        // Construct universe
        universes[parent_idx] = std::make_shared<Universe>(std::move(params));
    }
    CELER_ENSURE(std::all_of(universes.begin(), universes.end(), is_true{}));
    return universes;
}

//---------------------------------------------------------------------------//
/*!
 * Construct metadata
 */
auto GeometryBuilder::build_metadata(const MapProtoUniv& proto_univ_ids,
                                     VecBuiltProto&      build_results) const
    -> VecMetadata
{
    VecMetadata md;
    md.reserve(proto_univ_ids.size());

    for (auto& built : build_results)
    {
        md.push_back(std::move(built.md));
    }
    CELER_ENSURE(std::all_of(md.begin(), md.end(), is_true{}));
    return md;
}

//---------------------------------------------------------------------------//
/*!
 * Construct matids
 */
auto GeometryBuilder::build_matids(const VecBuiltProto& build_results,
                                   const VecUniverse&   universes) const
    -> VecMatid
{
    CELER_EXPECT(build_results.size() == universes.size());
    CELER_EXPECT(std::all_of(universes.begin(), universes.end(), is_true{}));

    VecMatid matids;

    // Concatenate matids from all universes
    for (auto univ_idx : range(build_results.size()))
    {
        auto local_num_volumes = universes[univ_idx]->num_volumes();
        auto local_vol_offset  = matids.size();
        // Resize to new maximum, filling with 'invalid'
        matids.resize(local_num_volumes + local_vol_offset,
                      geometria::invalid_matid());
        for (auto cell_matid : build_results[univ_idx].matids)
        {
            CELER_ASSERT(cell_matid.first.unchecked_get() < local_num_volumes);
            matids[cell_matid.first.get() + local_vol_offset]
                = cell_matid.second;
        }
    }
    return matids;
}

//---------------------------------------------------------------------------//
/*!
 * Construct boundary conditions
 */
auto GeometryBuilder::build_boundaries(const UniverseMetadata& top_md_base,
                                       MapShapeFaceBoundary boundaries) const
    -> MapSurfaceBoundary
{
    auto* top_md = dynamic_cast<const UnitMetadata*>(&top_md_base);
    Insist(top_md, "Top-level universe must be a Unit type");

    MapSurfaceBoundary result;

    // Loop over all surface IDs and each associated face
    for (auto surf_idx : range(top_md->num_surfaces()))
    {
        for (const ShapeFace& shape_face :
             top_md->surface_md(SurfaceId{surf_idx}))
        {
            auto iter = boundaries.find(shape_face);
            if (iter == boundaries.end())
            {
                // No user-specified boundary
                continue;
            }

            // Add boundary condition to the constructed list
            auto iter_insertion
                = result.insert({SurfaceId{surf_idx}, iter->second});

            // TODO: this could happen if the outer world has two shapes
            // that share a common surface (e.g. cylinder plus cuboid
            // sharing a +z face). We could check that the conditions are
            // consistent and warn, and only error if inconsistent.
            Insist(iter_insertion.second,
                   "Boundary " << to_string(shape_face)
                               << " is the same surface as "
                                  "an already-added boundary");

            // Remove boundary from user-provided list
            boundaries.erase(iter);
        }
    }

    // Check that the user-provided boundary map is now empty
    CELER_VALIDATE(
        boundaries.empty(),
        << "Some requested boundary surfaces do not exist in the global "
           "unit '"
        << top_md->metadata().name() << "':"
        << join(boundaries.begin(),
                boundaries.end(),
                ", ",
                [](const celeritas::GeometryBuilder::MapShapeFaceBoundary::value_type&
                       face_bc) { return to_string(face_bc.first); }));
    return result;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
