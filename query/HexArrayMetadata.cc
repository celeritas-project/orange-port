//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/HexArrayMetadata.cc
 * \brief HexArrayMetadata class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "HexArrayMetadata.hh"

#include <string>

#include "base/Casts.hh"
#include "base/Face.hh"
#include "base/Join.hh"
#include "base/Range.hh"
#include "base/StringFunctions.hh"
#include "orange/track/HexArrayTracker.hh"

namespace celeritas
{
namespace
{
//---------------------------------------------------------------------------//
/*!
 * Convert hex orientation to equivalent C string
 */
const char* to_cstring(HexArrayTracker::Orientation type)
{
#define ORANGE_ENUM_CASE(VAL)               \
    case HexArrayTracker::Orientation::VAL: \
        return #VAL

    switch (type)
    {
        ORANGE_ENUM_CASE(pointy_top);
        ORANGE_ENUM_CASE(flat_top);
    }

#undef ORANGE_ENUM_CASE
    CELER_ASSERT_UNREACHABLE();
}
} // namespace

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from cell/surface metadata
 */
HexArrayMetadata::HexArrayMetadata(Params params)
    : md_(std::move(params.unit))
    , cell_indexer_(params.dims)
    , bbox_(params.bbox)
{
    CELER_EXPECT(md_);
    CELER_EXPECT(bbox_);
}

//---------------------------------------------------------------------------//
//! Virtual destructor
HexArrayMetadata::~HexArrayMetadata() = default;

//---------------------------------------------------------------------------//
/*!
 * Get metadata about this universe
 */
const ObjectMetadata& HexArrayMetadata::metadata() const
{
    return md_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Local bounding box (extents) of this universe
 */
auto HexArrayMetadata::bbox() const -> BoundingBox
{
    return bbox_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Number of local surfaces
 */
auto HexArrayMetadata::num_surfaces() const -> surface_int
{
    return 8;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Number of local cells
 */
auto HexArrayMetadata::num_volumes() const -> volume_int
{
    return cell_indexer_.size();
}

//---------------------------------------------------------------------------//
/*!
 * Convert a local surface ID into a user-facing surface label
 */
std::string HexArrayMetadata::id_to_label(SurfaceId id) const
{
    CELER_EXPECT(id < this->num_surfaces());
    return to_string(HexprismFace<int>(id.get()));
}

//---------------------------------------------------------------------------//
/*!
 * Convert a local cell ID into a user-facing cell label
 */
std::string HexArrayMetadata::id_to_label(VolumeId id) const
{
    CELER_EXPECT(id < this->num_volumes());
    auto uvz = cell_indexer_.index(id.unchecked_get());

    std::ostringstream os;
    os << '{' << uvz[0] << ',' << uvz[1] << ',' << uvz[2] << '}';
    return os.str();
}

//---------------------------------------------------------------------------//
/*!
 * Get the volume of a cell
 */
real_type HexArrayMetadata::volume(VolumeId id) const
{
    CELER_EXPECT(id < cell_indexer_.size());
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Set the volume of a cell
 */
void HexArrayMetadata::set_volume(VolumeId id, real_type vol)
{
    CELER_EXPECT(id < this->num_volumes());
    CELER_EXPECT(vol >= 0);

    NotImplemented("Setting hexagonal array cell volumes");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Describe using data from the corresponding tracker
 */
void HexArrayMetadata::describe(std::ostream&  os,
                                const Tracker& base_tracker) const
{
    const auto& tracker = smart_cast<const HexArrayTracker&>(base_tracker);
    const auto& dims    = tracker.dims();
    const auto& z_edges = tracker.z_edges();

    os << ":Type: hexagonal array\n"
       << ":Orientation: " << to_cstring(tracker.orientation()) << '\n'
       << ":Apothem: " << tracker.apothem() << '\n'
       << ":Dimensions: [" << join(dims.begin(), dims.end(), ", ") << "]\n"
       << ":z: [" << join(z_edges.begin(), z_edges.end(), ", ") << "]\n\n";
}

//---------------------------------------------------------------------------//
} // namespace celeritas
