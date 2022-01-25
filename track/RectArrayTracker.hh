//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/RectArrayTracker.hh
 * \brief RectArrayTracker class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Tracker.hh"

#include "base/RegularGrid.hh"
#include "Definitions.hh"

namespace nemesis
{
template<class T>
struct CuboidFaceTraits;
}

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Track a particle on a cartesian brick mesh.
 *
 * TODO: add ArrayTracker hierarchy helper class, since many of the "special"
 * methods and types are identical between rect/hex/dode.
 */
class RectArrayTracker final : public Tracker
{
  public:
    //@{
    //! Public type aliases
    using Grid_t    = RegularGrid;
    using DimVector = Grid_t::DimVector;
    using Face      = ArrayFace<CuboidFaceTraits>;
    //@}

  public:
    //// CONSTRUCTION ////

    // Construct from cartesian grid
    explicit RectArrayTracker(Grid_t grid);

    //// TRACKING ////

    // Find the local cell and possibly surface ID.
    Initialization initialize(LocalState state) const final;

    // Calculate distance-to-intercept for the next surface
    Intersection intersect(LocalState state) const final;

    // Calculate normal on the current surface
    Real3 normal(LocalState state) const final;

    //// ACCESSORS ////

    //! Number of cells
    size_type num_volumes() const final { return grid_.num_volumes(); }

    //! Number of surfaces
    size_type num_surfaces() const final { return 6; }

    //! Access cell grid
    const Grid_t& grid() const { return grid_; }

    //! Cell id for given coordinate
    VolumeId volume_id(const DimVector& c) const { return coord_to_cell(c); }

  private:
    //// DATA ////

    Grid_t grid_;

    //// IMPLEMENTATION ////

    DimVector      cell_to_coords(VolumeId cell) const;
    VolumeId       coord_to_cell(const DimVector& ijk) const;
    Initialization initialize_interior(LocalState) const;
    Initialization cross_surface(DimVector coords, Face face) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
