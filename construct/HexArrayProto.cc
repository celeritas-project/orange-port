//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/HexArrayProto.cc
 * \brief HexArrayProto class definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "HexArrayProto.hh"

#include "base/Casts.hh"
#include "base/Assert.hh"
#include "base/SoftEquivalence.hh"
#include "base/Constants.hh"
#include "base/Definitions.hh"
#include "base/Face.hh"
#include "base/FixedViewArray.hh"
#include "base/Future.hh"
#include "base/Range.hh"
#include "orange/BoundingBox.hh"
#include "orange/Fuzziness.hh"
#include "orange/track/Definitions.hh"
#include "orange/track/HexArrayTracker.hh"
#include "orange/query/HexArrayMetadata.hh"
#include "PrismShape.hh"
#include "CuboidShape.hh"
#include "PlacedShape.hh"

using Axis::x;
using Axis::y;
using Axis::z;
using HexFace        = celeritas::HexArrayTracker::Face;
using HexOrientation = celeritas::HexArrayProto::Orientation;

namespace celeritas
{
namespace
{
//---------------------------------------------------------------------------//

constexpr int       U              = static_cast<int>(HexFace::Axes::u);
constexpr int       V              = static_cast<int>(HexFace::Axes::v);
constexpr real_type inv_sqrt_three = 0.5773502691896258;

//---------------------------------------------------------------------------//

struct HexEntryMd
{
    const ObjectMetadata&           parent_md;
    const HexArrayProto::DimVector& uvw;
    const ObjectMetadata&           daughter_md;
};

std::ostream& operator<<(std::ostream& os, const HexEntryMd& v)
{
    os << "hex array " << v.parent_md << " element {" << v.uvw << "}, "
       << v.daughter_md;
    return os;
}

//---------------------------------------------------------------------------//
/*!
 * Get orientation of a hex based on a prism shape.
 */
HexOrientation get_orientation(const PrismShape& prism)
{
    HexOrientation orientation = HexOrientation::automatic;
    real_type      rotation    = prism.rotation();
    CELER_ASSERT(rotation >= 0 && rotation < 1);

    if (soft_equiv(rotation, 0.0))
    {
        orientation = HexOrientation::flat_top;
    }
    else if (soft_equiv(rotation, half))
    {
        orientation = HexOrientation::pointy_top;
    }
    return orientation;
}

//---------------------------------------------------------------------------//
/*!
 * Convert hex orientation to equivalent C string
 */
const char* to_cstring(HexOrientation type)
{
#define ORANGE_ENUM_CASE(VAL) \
    case HexOrientation::VAL: \
        return #VAL

    switch (type)
    {
        ORANGE_ENUM_CASE(pointy_top);
        ORANGE_ENUM_CASE(flat_top);
        ORANGE_ENUM_CASE(automatic);
    }

#undef ORANGE_ENUM_CASE
    CELER_ASSERT_UNREACHABLE();
}

//---------------------------------------------------------------------------//
} // namespace

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
HexArrayProto::HexArrayProto(Params params)
    : md_(std::move(params.md))
    , units_(std::move(params.units))
    , layout_(params.layout)
    , orientation_(params.orientation)
{
    Insist(md_, "Missing array metadata");
    Insist(!units_.empty(),
           "No units are defined for array '" << md_.name() << "'");

    VecDbl heights(units_.dims()[Z], 0.0);
    apothem_ = 0;

    const auto& indexer = units_.indexer();
    for (auto user_volume_idx : range(indexer.size()))
    {
        DimVector           uvz     = indexer.index(user_volume_idx);
        const SPConstProto& sp_unit = units_[uvz];
        CELER_VALIDATE(
            sp_unit, << "In " << md_ << ": unit entry " << uvz << " is empty");
        const HexEntryMd err_where{md_, uvz, sp_unit->metadata()};

        // Check boundary entry's count, sense, and shape
        const RegionVec& boundaries = sp_unit->interior();
        Validate(boundaries.size() == 1
                     && boundaries.front().first == Sense::inside,
                 "In " << err_where << ": boundary must be a single hexagon");

        // Check that no extra rotation is present
        const PlacedShape& placed = *boundaries.front().second;
        CELER_VALIDATE(!placed.transform().has_rotation(),
                       << "In " << err_where
                       << ": cannot have transform rotation");

        // Boundary shape
        const auto* prism
            = dynamic_cast<const PrismShape*>(placed.shape().get());
        CELER_VALIDATE(prism,
                       << "In " << err_where << ": invalid bounding shape "
                       << placed.shape()->type());

        // Prism must have 6 side faces
        CELER_VALIDATE(prism->num_sides() == 6,
                       << "In " << err_where << ": invalid bounding shape has "
                       << prism->num_sides() << " sides");

        // Get prism properties and assign if unset
        auto      orientation = get_orientation(*prism);
        real_type height      = prism->upper() - prism->lower();
        real_type apothem     = prism->apothem();

        if (orientation_ == Orientation::automatic)
        {
            orientation_ = orientation;
        }
        if (heights[uvz[Z]] == 0)
        {
            heights[uvz[Z]] = height;
        }
        if (apothem_ == 0)
        {
            apothem_ = apothem;
        }

        CELER_VALIDATE(orientation != Orientation::automatic,
                       << "In " << err_where
                       << ": invalid prism rotation, must be 0 or half but is "
                       << prism->rotation());
        CELER_VALIDATE(orientation == orientation_,
                       << "In " << err_where
                       << ": inconsistent hex orientation, must be "
                       << to_cstring(orientation_));
        CELER_VALIDATE(apothem == apothem_,
                       << "In " << err_where
                       << ": inconsistent hex size, expected inner radius "
                          "(apothem) "
                       << apothem_ << " but got " << apothem);
        CELER_VALIDATE(height == heights[uvz[Z]],
                       << "In " << err_where
                       << ": inconsistent hex height, expected "
                       << heights[uvz[Z]] << " but got " << height);
    }

    // Integrate heights to become edges
    real_type start = 0.0;
    for (auto& pos : heights)
    {
        CELER_ASSERT(pos > 0);
        real_type cur_height = pos;
        pos                  = start;
        start += cur_height;
    }
    heights.push_back(start);
    z_edges_ = std::move(heights);

    // Axis-aligned (i) and off-axis (j) axes: U/X or V/Y
    i_ = (orientation_ == Orientation::pointy_top) ? U : V;
    j_ = (orientation_ == Orientation::pointy_top) ? V : U;

    // Hex array tracker expands rectangular layout on a rhomboidal
    // layout. Get the change in size from user -> tracker dimensions.
    DimVector tracker_dims = units_.dims();
    if (layout_ == Layout::rectangular)
    {
        i_offset_ = (tracker_dims[j_] - 1) / 2;
    }
    tracker_dims[i_] += i_offset_;
    tracker_indexer_ = ProtoVec3::Indexer_t(tracker_dims);

    // Calculate bounding box enclosing every cell
    Real3 lower = this->calc_origin({0, 0, 0});
    Real3 upper = this->calc_origin(indexer.dims() - DimVector{1, 1, 1});
    const real_type radius = 2 * inv_sqrt_three * apothem_;

    lower[i_] -= apothem_;
    lower[j_] -= radius;
    upper[i_] += apothem_;
    upper[j_] += radius;
    upper[Z] = z_edges_.back();
    if (layout_ == Layout::rectangular && units_.dims()[j_] % 2)
    {
        // The even hex sticks out more than the "last" hex along the
        // axis-aligned direction
        upper[i_] += apothem_;
    }

    CELER_ASSERT(lower.all_lt(upper));
    bbox_ = {lower, upper};

    CELER_ENSURE(z_edges_.size() > 1);
    CELER_ENSURE(z_edges_.front() == 0.0);
    CELER_ENSURE(apothem_ > 0.);
    CELER_ENSURE(orientation_ != Orientation::automatic);
    CELER_ENSURE(bbox_);
}

//---------------------------------------------------------------------------//
//! Default destructor
HexArrayProto::~HexArrayProto() = default;

//---------------------------------------------------------------------------//
/*!
 * Transformation to move the 'low' array corner to the specified point.
 *
 * This particular function is only used for the GG-style XML builder, not for
 * KENO.  The low array corner defined by GG (and used by ORANGE for
 * compatibility) is the *outer* bounding box for the hex array.
 *
 * The function returns the daughter-to-parent transform that will place the
 * outer lower-left external bounding corner at the position \c pos in the
 * parent unit.
 */
auto HexArrayProto::calc_placement(SpanConstReal3 pos) const -> Transform
{
    Real3 translation = -bbox_.lower();
    translation += make_vector(pos);
    return Transform(translation);
}

//---------------------------------------------------------------------------//
/*!
 * Transformation to move the center of the specified unit.
 *
 * The transform is the daughter-to-parent transform that will place the
 * *center* of array unit \c uvz at the position \c pos in the parent unit.
 *
 * This is not the same as the *origin* of the array unit.
 */
auto HexArrayProto::calc_placement(DimVector uvz, SpanConstReal3 pos) const
    -> Transform
{
    CELER_EXPECT(uvz.all_lt(units_.dims()));
    CELER_EXPECT(units_[uvz]);

    Real3 translation(-this->calc_translation(uvz));

    // Array coordinates => enclosing parent unit
    translation += make_vector(pos);

    return translation;
}

//---------------------------------------------------------------------------//
/*!
 * Construct a hex array from this proto.
 *
 * The construction method is tied closely to the implementation of the hex
 * array proto, specifically with cell ordering. The hex array always works in
 * logically rectangular (rhomboidal) layout, so we convert to that.
 */
auto HexArrayProto::build(BuildArgs args) const -> BuildResult
{
    CELER_EXPECT(args.implicit_boundary);

    using TrackerOrientation = HexArrayTracker::Orientation;

    // {u,v,z} dimensions of user-provided universes
    const auto& user_dims = units_.dims();

    // Construct the transformed daughters
    BuildResult result;

    // Fill daughters
    for (auto u : range(user_dims[U]))
    {
        for (auto v : range(user_dims[V]))
        {
            for (auto z : range(user_dims[Z]))
            {
                // Low corner of array cell => array coordinates
                DimVector uvz         = {u, v, z};
                Real3     translation = this->calc_translation(uvz);

                result.daughters.push_back(
                    {VolumeId(
                         tracker_indexer_.index(this->user_to_tracker(uvz))),
                     {units_[uvz], Transform(translation)}});
                CELER_ASSERT(result.daughters.back().first
                             < tracker_indexer_.size());
            }
        }
    }

    // Create tracker
    result.tracker = make_unique<HexArrayTracker>(
        apothem_,
        HexArrayTracker::PlaneDimVector(tracker_indexer_.dims()[U],
                                        tracker_indexer_.dims()[V]),
        z_edges_,
        orientation_ == Orientation::flat_top ? TrackerOrientation::flat_top
                                              : TrackerOrientation::pointy_top);

    // Create metadata
    HexArrayMetadata::Params md_params;
    md_params.unit = md_;
    md_params.dims = tracker_indexer_.dims();
    md_params.bbox = bbox_;
    result.md      = std::make_shared<HexArrayMetadata>(std::move(md_params));

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Get the boundary definition for defining a hole in a higher level
 *
 * TODO: calculate the largest hex shape (with similar flat/pointy to the array
 * definition) that fits entirely within the defined array cells.
 */
auto HexArrayProto::interior() const -> const RegionVec&
{
    NotImplemented("Implicitly creating a boundary from a hex array");
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the array cell ID of a user-given coordinate.
 */
auto HexArrayProto::user_to_tracker(const DimVector& user_uvz) const
    -> DimVector
{
    DimVector tracker_uvz = user_uvz;
    if (layout_ == Layout::rectangular)
    {
        // Transform from user to internal coordinates
        CELER_ASSERT(tracker_uvz[i_] + i_offset_ >= tracker_uvz[j_] / 2);
        tracker_uvz[i_] += i_offset_ - tracker_uvz[j_] / 2;
    }
    return tracker_uvz;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the base of a hex.
 *
 * The origin should always be at the bottom center of the lower-left *tracker*
 * hex (i.e. expanded to a rhomboid layout).
 */
Real3 HexArrayProto::calc_origin(const DimVector& uvz) const
{
    CELER_EXPECT(uvz.all_lt(units_.dims()));

    DimVector tracker_uvz = this->user_to_tracker(uvz);

    // Direction along normal (i) and off-normal (j) axes between hex centers
    const real_type wi = 2 * apothem_;
    const real_type wj = 3 * inv_sqrt_three * apothem_; // 1.5 * radius
    Real3           origin;
    origin[i_] = wi * tracker_uvz[i_] + apothem_ * tracker_uvz[j_];
    origin[j_] = wj * tracker_uvz[j_];
    origin[Z]  = z_edges_[tracker_uvz[Z]];

    return origin;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the transform from parent to the daughter unit of the hex.
 */
Real3 HexArrayProto::calc_translation(const DimVector& uvz) const
{
    CELER_EXPECT(uvz.all_lt(units_.dims()));

    using smart_cast;

    // Within-array origin of the cell: array origin is center bottom of
    // (0,0,0)
    Real3 trans = this->calc_origin(uvz);

    // Get "placed shape" and apply its translation
    const Proto::RegionVec& interior = units_[uvz]->interior();
    CELER_ASSERT(interior.size() == 1);
    CELER_ASSERT(interior.front().first == inside);
    const PlacedShape& ps = *interior.front().second;
    trans -= ps.transform().translation();

    // Get hex shape and subtract its translation
    const PrismShape& hex = smart_cast<const PrismShape&>(*ps.shape());
    trans[Z] -= hex.lower();

    return trans;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
