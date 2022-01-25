//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/UnitMetadata.hh
 * \brief UnitMetadata class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "UniverseMetadata.hh"

#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "Nemesis/containers/Span.hh"
#include "base/FlatTable.hh"
#include "orange/BoundingBox.hh"
#include "orange/construct/PlacedShape.hh"
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Map unit tracker cells and surfaces to construction metadata
 *
 * The return type for a surface is a view to pairs of shapes and surface
 * extensions.. If one surface is shared by multiple shapes, that surface ID
 * will correspond to multiple shape/extension pairs. So a plane might be a
 * vector that looks like:
 * \code
 * {{left_box, "px"}, {right_box, "mx"}}
 * \endcode
 *
 * to do:
 * - Don't reconstruct cells and displaying low-level logic (though we
 *   may still want to an option include that for testing?)
 * - Print mixed bottom-up CSG tree: surfaces, shapes, cells
 */
class UnitMetadata final : public UniverseMetadata
{
  public:
    //@{
    //! Types
    using SPConstShape       = std::shared_ptr<const PlacedShape>;
    using ShapeFace          = std::pair<SPConstShape, std::string>;
    using SpanConstShapeFace = span<const ShapeFace>;
    using TableShapeFace     = FlatTable<ShapeFace>;
    using VecMetadata        = std::vector<ObjectMetadata>;
    using VecDbl             = std::vector<real_type>;
    //@}

    struct Params
    {
        ObjectMetadata unit;
        TableShapeFace surfaces;
        VecMetadata    cells;
        BoundingBox    bbox;
        VecDbl         volumes;
        bool           is_simple{false};
    };

  public:
    // Construct from cell/surface metadata
    UnitMetadata(Params params);

    // Virtual destructor
    virtual ~UnitMetadata();

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

    //// ACCESSORS ////

    //! Whether meaninigful volumes are all in the same Z order
    bool is_simple() const { return md_.is_simple; }

    // Get detailed metadata for a surface
    inline SpanConstShapeFace surface_md(SurfaceId) const;

    // Get detailed metadata for a cell
    inline const ObjectMetadata& vol_md(VolumeId) const;

  private:
    // Stored metadata
    Params md_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "UnitMetadata.i.hh"
//---------------------------------------------------------------------------//
