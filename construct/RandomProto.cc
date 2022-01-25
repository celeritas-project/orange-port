//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RandomProto.cc
 * \brief RandomProto class definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RandomProto.hh"

#include <algorithm>
#include <cstdint>
#include <memory>
#include "UnitProto.hh"
#include "orange/construct/UnitBuilder.hh"
#include "orange/track/MaskedUnitTracker.hh"
#include "Nemesis/database/Std_DB.hh"
#include "Nemesis/comm/Logger.hh"
#include "base/GeometryUtils.hh"
#include "base/Range.hh"
#include "base/Future.hh"
#include "Nemesis/sprng/Random.hh"
#include "Nemesis/serialize/Hasher.hh"
#include "orange/Definitions.hh"
#include "../construct/Shape.hh"
#include "../construct/SphereShape.hh"
#include "../packing/BuildPacker.hh"

#include <limits>

using RNGEngine  = random::RNGEngine;
using RNGControl = random::RNGControl;

namespace
{
//---------------------------------------------------------------------------//
/*!
 * Build an RNG based on a name and user-provided seed.
 *
 * This ensures that reusing a "seed" for multiple random protos will not make
 * them correlated.
 */
RNGEngine build_rng(const std::string& name, int seed)
{
    CELER_EXPECT(RNGControl::initialized());

    static const int max_int   = std::numeric_limits<int>::max();
    int              stream_id = 0;
    {
        // Hash the name, using 'xor' to reduce the 256-bit hash to 32-bit
        Hasher h;
        h << name;

        std::uint32_t combined = 0;
        for (std::uint32_t val : h.hash())
        {
            combined ^= val;
        }
        stream_id = (combined & static_cast<std::uint32_t>(max_int));
    }

    // XOR against user-given seed
    stream_id ^= seed;

    // Ensure the result is within the available number of streams
    const RNGControl& control = RNGControl::instance();
    if (control.max_streams() < max_int - 1)
    {
        stream_id = stream_id % (control.max_streams() + 1);
    }
    CELER_ASSERT(stream_id < control.max_streams());

    return control.engine(stream_id);
}

//---------------------------------------------------------------------------//
} // namespace

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Constructor
 */
RandomProto::RandomProto(RandomProto::Params params) : data_(std::move(params))
{
    Insist(data_.md, "Missing unit metadata");
    Insist(!data_.interior.empty(),
           "Missing interior definition for unit '" << data_.md.name() << "'");
    Insist(!data_.particles.empty(),
           "Missing particles in unit '" << data_.md.name() << "'");
    Insist(std::all_of(data_.particles.begin(),
                       data_.particles.end(),
                       [](const Particle& p) { return static_cast<bool>(p); }),
           "Incomplete particle definitions");
    Insist(data_.interior.size() == 1 && data_.interior.front().second
               && data_.interior.front().first == Sense::inside,
           "Currently random proto must have exactly one shape (with 'inside' "
           "sense)");
    const auto& boundary_shape = data_.interior.front().second->shape();
    CELER_VALIDATE(boundary_shape->inradius() > 0,
                   << "Invalid shape type '" << boundary_shape->type()
                   << "' for "
                      "'random' unit '"
                   << data_.md.name() << "'");

    for (const auto& p : data_.particles)
    {
        const auto& pint = p.proto->interior();
        CELER_VALIDATE(pint.size() == 1 && pint.front().first == Sense::inside,
                       << "Daughter '" << p.proto->metadata().name()
                       << "' cannot have a boundary with multiple "
                          "regions in random unit '"
                       << data_.md.name() << "'");
    }
}

//---------------------------------------------------------------------------//
//! Default destructor
RandomProto::~RandomProto() = default;

//---------------------------------------------------------------------------//
/*!
 * Build the universe.
 */
auto RandomProto::build(BuildArgs args) const -> BuildResult
{
    // Construct particle input from our stored protos
    std::vector<ParticleInput> pinput(data_.particles.size());
    for (auto i : range(pinput.size()))
    {
        const auto& proto = *data_.particles[i].proto;

        // Access the boundary's raw shape
        const Shape& shape = *proto.interior().front().second->shape();
        const auto*  sph   = dynamic_cast<const SphereShape*>(&shape);
        CELER_VALIDATE(sph,
                       << "Cannot add non-spheres to the 'random' universe '"
                       << data_.md.name() << "': embedded universe '"
                       << proto.metadata().name() << "' has bounding shape '"
                       << shape.type() << "'");

        // Construct particle input
        ParticleInput& pi  = pinput[i];
        pi.volume_fraction = data_.particles[i].volume_fraction;
        pi.label           = i;
        pi.radius          = sph->radius();
    }

    // Construct particle type
    const auto& boundary_shape = *data_.interior.front().second->shape();
    auto        pt             = std::make_shared<ParticleTypes>(
        std::move(pinput), boundary_shape.volume(), boundary_shape.inradius());

    // Construct RNG if it hasn't yet been done
    if (!RNGControl::initialized())
    {
        RNGControl::initialize(20171213);
    }

    // Construct packer
    auto db = std::make_shared<database::Std_DB>("packer_db");
    db->set<int>("failure_batch_size", data_.options.failure_batch_size);
    db->set<real_type>("failure_tolerance", data_.options.failure_tolerance);

    std::unique_ptr<PackerBase> packer
        = build_packer(boundary_shape,
                       build_rng(data_.md.name(), data_.options.seed),
                       std::move(db));
    CELER_ASSERT(packer);

    // Pack the particles
    log() << "Instantiating random universe '" << data_.md.name() << "'";
    ParticleCollection particles = packer->pack_random_rejection(pt);

    // Construct a unit with particles as holes
    BuildResult result;
    UnitBuilder build_unit;

    build_unit.exterior(data_.interior, ZOrder::implicit_exterior, data_.md);
    real_type fill_media_volume = boundary_shape.volume();

    // Loop over particles, constructing attributes of the different particle
    for (auto ptype_idx : range(particles.num_types()))
    {
        // Note that particle types may not be in the same order as protos_,
        // so use the label to get the correct one. Get the universe
        // corresponding to this
        CELER_ASSERT(pt->label(ptype_idx) < data_.particles.size());

        // Attributes of the particles of this type
        auto      pt_centers = particles.centers(ptype_idx);
        real_type pt_radius  = particles.ptypes()->radius(ptype_idx);

        // Subtract the total particle volume from the interior volume
        fill_media_volume -= sphere_volume(pt_radius) * pt_centers.size();

        // Place particles via holes
        const Particle& p = data_.particles[pt->label(ptype_idx)];
        for (const Real3& pos : make_span(pt_centers))
        {
            UnitProto::Hole particle_hole;
            particle_hole.proto     = p.proto;
            particle_hole.transform = Transform(pos);
            particle_hole.md        = p.md;
            SPConstShape hole_boundary_shape
                = UnitProto::make_hole_shape(particle_hole);

            // Add region
            auto volume_id = build_unit.region(
                {{neg, hole_boundary_shape}}, ZOrder::hole, particle_hole.md);

            // Add daughter universe
            result.daughters.push_back(
                {volume_id,
                 Daughter{particle_hole.proto, particle_hole.transform}});
        }
    }

    // Build the fill media
    {
        auto volume_id = build_unit.region(
            data_.interior, ZOrder::media, data_.md, fill_media_volume);
        result.matids.push_back({volume_id, data_.fill_matid});
    }

    // Construct tracker: use masked tracker because we're using multiple
    // zorders (holes overried background medium)
    auto unit_components = build_unit(data_.md);
    result.tracker       = make_unique<MaskedUnitTracker>(
        std::move(unit_components.surfaces), unit_components.regions);
    result.md = std::move(unit_components.md);
    return result;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
