//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/RectArrayMetadata.cc
 * \brief RectArrayMetadata class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RectArrayMetadata.hh"

#include "base/Casts.hh"
#include "base/Face.hh"
#include "base/Join.hh"
#include "base/Range.hh"
#include "base/StringFunctions.hh"
#include "orange/track/RectArrayTracker.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct from cell/surface metadata
 */
RectArrayMetadata::RectArrayMetadata(Params params)
    : md_(std::move(params.unit))
    , cell_indexer_({params.dims[2], params.dims[1], params.dims[0]})
    , bbox_(params.bbox)
{
}

//---------------------------------------------------------------------------//
//! Virtual destructor
RectArrayMetadata::~RectArrayMetadata() = default;

//---------------------------------------------------------------------------//
/*!
 * Get metadata about this universe
 */
const ObjectMetadata& RectArrayMetadata::metadata() const
{
    return md_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Local bounding box (extents) of this universe
 */
auto RectArrayMetadata::bbox() const -> BoundingBox
{
    return bbox_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Number of local surfaces
 */
auto RectArrayMetadata::num_surfaces() const -> surface_int
{
    return 6;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Number of local cells
 */
auto RectArrayMetadata::num_volumes() const -> volume_int
{
    return cell_indexer_.size();
}

//---------------------------------------------------------------------------//
/*!
 * Convert a local surface ID into a user-facing surface label.
 */
std::string RectArrayMetadata::id_to_label(SurfaceId id) const
{
    CELER_EXPECT(id < this->num_surfaces());
    return to_string(CuboidFace<int>(id.get()));
}

//---------------------------------------------------------------------------//
/*!
 * Convert a local cell ID into a user-facing cell label
 */
std::string RectArrayMetadata::id_to_label(VolumeId id) const
{
    CELER_EXPECT(id < cell_indexer_.size());

    auto kji = cell_indexer_.index(id.get());

    std::ostringstream os;
    os << '{' << kji[2] << ',' << kji[1] << ',' << kji[0] << '}';
    return os.str();
}

//---------------------------------------------------------------------------//
/*!
 * Get the volume of a cell
 */
real_type RectArrayMetadata::volume(VolumeId id) const
{
    CELER_EXPECT(id < cell_indexer_.size());
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Set the volume of a cell
 */
void RectArrayMetadata::set_volume(VolumeId id, real_type vol)
{
    CELER_EXPECT(id < cell_indexer_.size());
    CELER_EXPECT(vol >= 0);

    NotImplemented("Setting rectangular array cell volumes");
}

//---------------------------------------------------------------------------//
/*!
 * Describe using data from the corresponding tracker
 */
void RectArrayMetadata::describe(std::ostream&  os,
                                 const Tracker& base_tracker) const
{
    const auto& tracker = smart_cast<const RectArrayTracker&>(base_tracker);

    auto dims = tracker.grid().num_cells_dims();
    os << ":Dimensions: [" << join(dims.begin(), dims.end(), ", ") << "]\n";
    for (auto ax : {Axis::x, Axis::y, Axis::z})
    {
        // Restore un-bumped boundaries to grid
        auto edges    = tracker.grid().edges(ax);
        edges.front() = bbox_.lower()[ax];
        edges.back()  = bbox_.upper()[ax];
        os << ':' << def::xyz_name(ax) << ": ["
           << join(edges.begin(), edges.end(), ", ") << "]\n";
    }
    os << '\n';
}

//---------------------------------------------------------------------------//
} // namespace celeritas
