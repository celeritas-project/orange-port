//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/DodeArrayMetadata.cc
 * \brief DodeArrayMetadata class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "DodeArrayMetadata.hh"

#include "base/Casts.hh"
#include "base/Face.hh"
#include "base/Join.hh"
#include "base/Range.hh"
#include "base/StringFunctions.hh"
#include "orange/track/DodeArrayTracker.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct from cell/surface metadata
 */
DodeArrayMetadata::DodeArrayMetadata(Params params)
    : md_(std::move(params.unit))
    , cell_indexer_(params.dims)
    , bbox_(params.bbox)
{
}

//---------------------------------------------------------------------------//
//! Virtual destructor
DodeArrayMetadata::~DodeArrayMetadata() = default;

//---------------------------------------------------------------------------//
/*!
 * Get metadata about this universe
 */
const ObjectMetadata& DodeArrayMetadata::metadata() const
{
    return md_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Local bounding box (extents) of this universe
 */
auto DodeArrayMetadata::bbox() const -> BoundingBox
{
    return bbox_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Number of local surfaces
 */
auto DodeArrayMetadata::num_surfaces() const -> surface_int
{
    return 12;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Number of local cells
 */
auto DodeArrayMetadata::num_volumes() const -> volume_int
{
    return cell_indexer_.size();
}

//---------------------------------------------------------------------------//
/*!
 * Convert a local surface ID into a user-facing surface label
 *
 * \todo Add better mapping to mesh surface indexer
 */
std::string DodeArrayMetadata::id_to_label(SurfaceId id) const
{
    return to_string(RhombicDodecahedronFace<int>(id.get()));
}

//---------------------------------------------------------------------------//
/*!
 * Convert a local cell ID into a user-facing cell label
 *
 * Note that unlike rect array which is indexed as KJI, dode array is indexed
 * naturally as IJK.
 */
std::string DodeArrayMetadata::id_to_label(VolumeId id) const
{
    CELER_EXPECT(id < cell_indexer_.size());

    auto ijk = cell_indexer_.index(id.get());

    std::ostringstream os;
    os << '{' << ijk[0] << ',' << ijk[1] << ',' << ijk[2] << '}';
    return os.str();
}

//---------------------------------------------------------------------------//
/*!
 * Get the volume of a cell
 */
real_type DodeArrayMetadata::volume(VolumeId id) const
{
    CELER_EXPECT(id < cell_indexer_.size());
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * Set the volume of a cell
 */
void DodeArrayMetadata::set_volume(VolumeId id, real_type vol)
{
    CELER_EXPECT(id < cell_indexer_.size());
    CELER_EXPECT(vol >= 0);

    NotImplemented("Setting dodecahedral array cell volumes");
}

//---------------------------------------------------------------------------//
/*!
 * Describe using data from the corresponding tracker
 */
void DodeArrayMetadata::describe(std::ostream&  os,
                                 const Tracker& base_tracker) const
{
    const auto& tracker = smart_cast<const DodeArrayTracker&>(base_tracker);
    const auto& dims    = tracker.dims();

    os << ":Type: dodecahedral array\n"
       << ":Apothem: " << tracker.apothem() << '\n'
       << ":Dimensions: [" << join(dims.begin(), dims.end(), ", ") << "]\n\n";
}

//---------------------------------------------------------------------------//
} // namespace celeritas
