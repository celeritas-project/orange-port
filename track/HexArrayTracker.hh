//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/HexArrayTracker.hh
 * \brief HexArrayTracker class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Tracker.hh"

#include "base/Definitions.hh"
#include "base/HexGrid.hh"

#include "orange/track/Definitions.hh"
#include "orange/surfaces/SurfaceContainer.hh"
#include "orange/BoundingBox.hh"

#include "Definitions.hh"

namespace nemesis
{
template<class T>
struct HexprismFaceTraits;
}

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Track a particle on a hexagonal grid.
 *
 * The internal tracking is always done on a rhomboidal array, which is
 * logically regular in UVZ (same number of cells in each row or column). The
 * origin of the array is the bottom of the center of the lower-left hex (u=0,
 * v=0, z=0), so the valid extents of the array will extend slightly negative
 */
class HexArrayTracker final : public Tracker
{
  public:
    //@{
    //! Public type aliases
    using volume_int     = VolumeId::size_type;
    using VecDbl         = std::vector<double>;
    using PlaneDimVector = Array<volume_int, 2>;
    using DimVector      = Array<volume_int, 3>;
    using PlaneVector    = Array<real_type, 2>;
    using Real3          = Array<real_type, 3>;
    using Indexer        = RegularIndexer<volume_int, 3>;
    using Face           = ArrayFace<HexprismFaceTraits>;
    //@}

    //! Orientation of the cells in the hex array.
    enum class Orientation
    {
        pointy_top, //!< U is along X
        flat_top    //!< V is along Y
    };

  public:
    //// CONSTRUCTION ////

    // Construct hex array tracker with the given data
    explicit HexArrayTracker(real_type      apothem,
                             PlaneDimVector uv_dims,
                             VecDbl         z_edges,
                             Orientation    orientation);

    //// TRACKING ////

    // Find the local cell and possibly surface ID.
    Initialization initialize(LocalState state) const final;

    // Calculate distance-to-intercept for the next surface
    Intersection intersect(LocalState state) const final;

    // Calculate normal on the current surface
    Real3 normal(LocalState state) const final;

    //// ACCESSORS ////

    //! Number of surfaces
    size_type num_surfaces() const final { return 8; }

    //! Number of cells
    size_type num_volumes() const final { return indexer_.size(); }

    //// TESTING HELPER FUNCTIONS ////

    //! Hex (u, v, z) dimensions
    const DimVector& dims() const { return indexer_.dims(); }

    // Cell id for given uvz coordinate
    VolumeId volume_id(const DimVector& c) const;

    //! Hex apothem (interior radius)
    real_type apothem() const { return apothem_; }

    // Circumradius (exterior radius)
    real_type circumradius() const;

    //! Flat or pointy tops
    Orientation orientation() const { return orientation_; }

    //! Z-planes of hex grid.
    const VecDbl& z_edges() const { return z_edges_; }

    // Find the logical index for the given point
    DimVector find(SpanConstReal3 point) const;

    // Return hex centroid (x,y,z) given (u,v,z) coordinates.
    Real3 centroid(const DimVector& coords) const;

  private:
    //// TYPES ////

    using SenseVector    = Array<Sense, 6>;
    using HexFaceNormals = Array<PlaneVector, 3>;

    //// DATA ////

    //! Inner radius of hex (apothem)
    real_type apothem_;

    //! Z grid
    VecDbl z_edges_;

    //! Orientation
    Orientation orientation_;

    //! Hex face normals for the first three faces
    HexFaceNormals uvw_normals_;

    //! Normal to the upward-pointing 3 faces divided by apothem length
    HexFaceNormals uvw_;

    //! (x,y) distances between hex centers along (u,v) directions
    Array<PlaneVector, 2> uv_span_;

    //! (U,V,Z) Indexer for hexagonal grid
    Indexer indexer_;

    //! Vector of six hex surfaces (not Z) in Face_t ordering
    SurfaceContainer surfaces_;

    // Vector of "inside" senses on each face
    SenseVector senses_;

    //// METHODS ////

    Initialization initialize_interior(LocalState) const;
    Initialization cross_surface(DimVector uvz, Face face) const;

    PlaneVector centroid(const PlaneDimVector& coords) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
