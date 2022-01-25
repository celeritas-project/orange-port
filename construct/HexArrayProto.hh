//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/HexArrayProto.hh
 * \brief HexArrayProto class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "ArrayProto.hh"

#include <vector>
#include "base/Definitions.hh"
#include "orange/BoundingBox.hh"
#include "orange/query/ObjectMetadata.hh"
#include "orange/track/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Define a hexagonal array of units.
 *
 * The resulting tracker will be in a rhombic configuration with the center
 * bottom of the (0,0,0) unit at the origin.
 */
class HexArrayProto final : public ArrayProto
{
  public:
    //! Logical offsetting of hex array.
    enum class Layout
    {
        rhomboidal, //!< rhomboidal, each row has the same length
        rectangular //!< rectangular, each row has the same length
    };

    //! Orientation of the cells in the hex array.
    enum class Orientation
    {
        pointy_top, //!< U is along X
        flat_top,   //!< V is along Y
        automatic   //!< auto-detect
    };

    //! Construction arguments
    struct Params
    {
        //! Orientation of hex array (pointy/flat/auto)
        Orientation orientation = Orientation::automatic;
        //! Logical offsetting of hex array (rhomb/rect)
        Layout layout = Layout::rhomboidal;
        //! {u,v,z} array entries from lower left
        ProtoVec3 units;
        //! Name and provenance
        ObjectMetadata md;
    };

  public:
    // Constructor
    explicit HexArrayProto(Params params);

    // Default destructor
    ~HexArrayProto();

    //! Get the proto metadata
    const ObjectMetadata& metadata() const final { return md_; }

    // Construct a boundary definition for defining a hole
    const RegionVec& interior() const final;

    // Daughter-to-parent transformation for the array lower-left corner
    Transform calc_placement(SpanConstReal3 pos) const final;

    // Transform to place the center of unit `index` in parent `pos`
    Transform calc_placement(DimVector index, SpanConstReal3 pos) const final;

    // Construct a universe from this proto
    BuildResult build(BuildArgs args) const final;

    //// HELPER FUNCTIONS ////

    // Convert from user to internal tracker cell ID
    DimVector user_to_tracker(const DimVector& uvz) const;

  private:
    //// TYPES ////

    using VecDbl = std::vector<double>;

    //// DATA ////

    ObjectMetadata md_;
    ProtoVec3      units_;
    real_type      apothem_;
    VecDbl         z_edges_;
    Layout         layout_;
    Orientation    orientation_;
    BoundingBox    bbox_;

    // Index translation and origin offsets
    ProtoVec3::Indexer_t tracker_indexer_;
    int                  i_        = -1; // U/X if pointy else V/Y
    int                  j_        = -1; // V/Y if pointy else U/X
    int                  i_offset_ = 0;  // Offset along i if rect

    //// METHODS ////

    // Validate the bounding shape on an array element
    void validate_cell_shape(const DimVector& uvz) const;

    // Add cell to array
    void calculate_array_parameters(const DimVector& uvz);

    // Calculate the origin of the given array coordinate
    Real3 calc_origin(const DimVector& uvz) const;

    // Calculate translation to the daughter element at uvz
    Real3 calc_translation(const DimVector& uvz) const;

    // Convert from user to internal tracker cell ID
    VolumeId volume_id(const DimVector& uvz) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas
