//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/DodeArrayProto.cc
 * \brief DodeArrayProto class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "DodeArrayProto.hh"

#include <memory>
#include <unordered_set>
#include "base/Casts.hh"
#include "base/SoftEquivalence.hh"
#include "base/Constants.hh"
#include "base/Definitions.hh"
#include "base/FixedViewArray.hh"
#include "base/Future.hh"
#include "base/Range.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "orange/track/DodeArrayTracker.hh"
#include "orange/query/DodeArrayMetadata.hh"
#include "RhombicDodecahedronShape.hh"
#include "CuboidShape.hh"
#include "PlacedShape.hh"

using make_unique;
using smart_cast;
using Axis::x;
using Axis::y;
using Axis::z;
using constants::sqrt_two;
using geometria::BoundingBox;

namespace celeritas
{
namespace
{
//---------------------------------------------------------------------------//
/*!
 * Obtain encompasing axis-aligned bounding box for dodecahedral array.
 */
BoundingBox
calc_bounding_box(real_type apothem, const DodeArrayProto::DimVector& dims)
{
    // This fully encompasses the dodecahedral array which means there are
    // undefined spaces on each face
    Real3 lower, upper;

    for (int ax : {X, Y})
    {
        lower[ax] = -apothem;
        upper[ax] = lower[ax] + apothem * 2 * dims[ax];
    }

    lower[Z]         = -sqrt_two * apothem;
    auto cell_height = 2.0 * sqrt_two * apothem;
    upper[Z]         = lower[Z] + cell_height * std::floor(dims[Z] / 2.0)
               + (dims[Z] % 2) * cell_height;

    CELER_ENSURE(lower < upper);
    return BoundingBox(lower, upper);
}
} // namespace

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Constructor
 */
DodeArrayProto::DodeArrayProto(Params params)
    : md_(std::move(params.md)), units_(std::move(params.units))
{
    Insist(md_, "Missing array metadata");
    Insist(!units_.empty(),
           "No units are defined for array '" << md_.name() << "'");

    std::unordered_set<const Proto*> checked_protos;
    apothem_ = -1;

    // Set apothem and validate all protos
    for (const SPConstProto& sp_unit : units_)
    {
        Insist(sp_unit, md_.name() << " has an empty unit");

        // See if we need to check the daughter unit
        auto iter_inserted = checked_protos.insert(sp_unit.get());
        if (!iter_inserted.second)
            continue;

        const auto& proto_name = sp_unit->metadata().name();

        // Check boundary entry's count and sens
        CELER_VALIDATE(sp_unit->interior().size() == 1,
                       << "Dodecahedral array " << md_.name() << " entry "
                       << proto_name
                       << " boundary must be a single dodecahedron");
        const auto& ss = sp_unit->interior().front();
        CELER_VALIDATE(ss.first == Sense::inside,
                       << "Dodecahedral array " << md_.name() << " entry "
                       << proto_name
                       << " boundary must be inside a single dodecahedron");

        // Check that the first external surface is a rhombic
        // dodecahedron
        const auto* rd = dynamic_cast<const RhombicDodecahedronShape*>(
            ss.second->shape().get());
        CELER_VALIDATE(rd,
                       << "Cannot fill an array element of '" << md_.name()
                       << "' with a non-dodecahedron shape: '" << proto_name
                       << "' has bounding shape '"
                       << ss.second->shape()->type() << "'");

        if (apothem_ <= 0)
        {
            // First unit: construct apothem
            apothem_ = rd->apothem();
        }

        // Check consistency with initial size (after initialization)
        CELER_VALIDATE(rd->apothem() == apothem_,
                       << "Dodecahedral array " << md_.name() << " entry "
                       << proto_name << " has inconsistent radii: "
                       << rd->apothem() << " vs " << apothem_);
    }
    CELER_ENSURE(apothem_ > 0);
}

//---------------------------------------------------------------------------//
//! Default destructor
DodeArrayProto::~DodeArrayProto() = default;

//---------------------------------------------------------------------------//
/*!
 * Transformation to move the 'low' array corner to the specified point.
 *
 * This is the daughter-to-parent transform: destination in enclosing unit
 * (function argument) destination; and array lower is always the origin.
 */
auto DodeArrayProto::calc_placement(SpanConstReal3 pos) const -> Transform
{
    Real3 translation = make_vector(pos);
    // Don't need to subtract array "origin" since it's always (0, 0, 0)
    return Transform(translation);
}

//---------------------------------------------------------------------------//
/*!
 * Transformation to move the center of the specified unit.
 *
 * The transform is the daughter-to-parent transform that will place the
 * *center* of array unit \c ijk at the position \c pos in the parent unit.
 *
 * This is not the same as the *origin* of the array unit.
 */
auto DodeArrayProto::calc_placement(DimVector ijk, SpanConstReal3 pos) const
    -> Transform
{
    CELER_EXPECT(ijk.all_lt(units_.dims()));
    CELER_EXPECT(units_[ijk]);

    Real3 translation(-calc_origin(ijk));

    // Array coordinates => enclosing parent unit
    translation += make_vector(pos);

    return translation;
}

Real3 DodeArrayProto::calc_origin(const DimVector& coords) const
{
    Real3 cell_origin(2.0 * apothem_ * coords[X],
                      2.0 * apothem_ * coords[Y],
                      apothem_ * coords[Z] * sqrt_two);

    // Correct X and Y for odd Z index array
    if (coords[Z] % 2 == 1)
    {
        cell_origin[X] += apothem_;
        cell_origin[Y] += apothem_;
    }
    return cell_origin;
}

//---------------------------------------------------------------------------//
/*!
 * Construct a universe from this proto
 */
auto DodeArrayProto::build(BuildArgs args) const -> BuildResult
{
    CELER_EXPECT(args.implicit_boundary);

    auto dims = units_.dims();

    // Construct the transformed daughters
    BuildResult result;
    result.daughters.reserve(units_.size());

    // Fill daughters
    for (auto k : range(dims[def::K]))
    {
        for (auto j : range(dims[def::J]))
        {
            for (auto i : range(dims[def::I]))
            {
                DimVector           coords  = {i, j, k};
                const SPConstProto& sp_unit = units_[coords];

                // Daughter => lower corner of array cell (0,0,0)
                Real3 translation = calc_origin(coords);

                VolumeId id(units_.indexer().index(coords));
                result.daughters.push_back(
                    {id, {sp_unit, Transform(translation)}});
            }
        }
    }

    // Create tracker
    result.tracker = make_unique<DodeArrayTracker>(apothem_, dims);

    // Create metadata
    DodeArrayMetadata::Params md_params;
    md_params.apothem = apothem_;
    md_params.unit    = md_;
    md_params.dims    = dims;
    md_params.bbox    = calc_bounding_box(apothem_, dims);
    result.md = std::make_shared<DodeArrayMetadata>(std::move(md_params));

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Get the boundary definition for defining a hole in a higher level
 */
auto DodeArrayProto::interior() const -> const RegionVec&
{
    NotImplemented("Implicitly creating a boundary from a dodecahedron array");
}

//---------------------------------------------------------------------------//
} // namespace celeritas
