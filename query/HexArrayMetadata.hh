//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/HexArrayMetadata.hh
 * \brief HexArrayMetadata class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "UniverseMetadata.hh"

#include "base/RegularIndexer.hh"
#include "orange/BoundingBox.hh"
#include "orange/Definitions.hh"
#include "orange/track/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Map hexagonal array cells and surfaces to names.
 */
//---------------------------------------------------------------------------//
class HexArrayMetadata final : public UniverseMetadata
{
  public:
    //@{
    //! Public type aliases
    using DimVector = Array<size_type, 3>;
    //@}

    struct Params
    {
        ObjectMetadata unit; //!< Unit metadata
        DimVector      dims; //!< {u, v, z} cells (rhomboid layout)
        BoundingBox    bbox; //!< Extents of mesh
    };

  public:
    // Constructor
    HexArrayMetadata(Params params);

    // Default destructor
    ~HexArrayMetadata();

    //// PUBLIC INTERFACE ////

    // Get metadata about this universe
    const ObjectMetadata& metadata() const final;

    // Local bounding box
    BoundingBox bbox() const final;

    // Number of local surfaces
    surface_int num_surfaces() const final;

    // Number of local cells
    volume_int num_volumes() const final;

    // Convert a local surface ID into a user-facing surface label
    std::string id_to_label(SurfaceId) const final;

    // Convert a local cell ID into a user-facing cell label
    std::string id_to_label(VolumeId) const final;

    // Get the volume of a cell
    real_type volume(VolumeId) const final;

    // Set the the volume of a cell
    void set_volume(VolumeId id, real_type volume) final;

    // Describe using data from the corresponding tracker
    void describe(std::ostream& os, const Tracker& tracker) const final;

  private:
    using CellIndexer_t = RegularIndexer<size_type, 3>;

    //// DATA ////

    ObjectMetadata md_;
    CellIndexer_t  cell_indexer_;
    BoundingBox    bbox_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
