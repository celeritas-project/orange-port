//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/RectArrayMetadata.hh
 * \brief RectArrayMetadata class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "UniverseMetadata.hh"

#include "base/RegularGrid.hh"
#include "orange/BoundingBox.hh"
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Map rectangular array cells and surfaces to names.
 *
 * \todo Currently cell volumes aren't implemented, but if we allow array cells
 * to be filled with a material then we'll need to take a grid rather than just
 * a set of dimensions.
 */
//---------------------------------------------------------------------------//
class RectArrayMetadata final : public UniverseMetadata
{
  public:
    //@{
    //! Public type aliases
    using DimVector = Array<size_type, 3>;
    //@}

    struct Params
    {
        ObjectMetadata unit; //!< Unit metadata
        DimVector      dims; //!< X, Y, Z cells
        BoundingBox    bbox; //!< Extents of mesh
    };

  public:
    // Constructor
    RectArrayMetadata(Params params);

    // Default destructor
    ~RectArrayMetadata();

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
