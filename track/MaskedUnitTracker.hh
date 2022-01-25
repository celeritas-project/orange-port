//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/MaskedUnitTracker.hh
 * \brief MaskedUnitTracker class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "UnitTracker.hh"

#include "detail/SurfaceGroupContainer.hh"
#include "detail/UnitTrackingStorage.hh"

namespace celeritas
{
struct UnitRegion;
//---------------------------------------------------------------------------//
/*!
 * Track a particle in a universe of overlapping cells
 *
 * The masked unit tracker is based on a set of surfaces and region definitions
 * based on those surfaces (usually derived from shapes).. Each region
 * definition is given a "Z" order that determines the precedence of the
 * region: whether it obscures a region of lower order while tracking. For
 * example, holes have a higher z order than media and thus take precedence.
 *
 * A kd-tree acceleration structure is built based on the bounding boxes
 * provided by each region.
 *
 * Cell IDs are local and correspond to the regions in the input vector. The
 * provided vector *must* be sorted in descending Z order, so that
 * higher-precedence regions are earlier in the list. This means the "exterior"
 * (implicit complement of the "interior" region definition vector) *MUST* be
 * the *first* region in the list!
 *
 * The tracker *never* modifies the state: it only takes local information
 * about the state and returns the change to be performed. The calling class is
 * responsible for updating spatial positions, changing cells, etc.
 *
 * The \c initialize routine is called when first initializing a particle
 * (looking for the enclosing cell) \em and when entering a new cell
 * after crossing a surface. In the former case, the local state will have
 * neither a cell nor surface ID; in the latter case, it will have \em both.
 * During initialization, all cells near the state's position are checked for
 * overlap. Overlap is defined as two cells with the same Z order precedence
 * sharing the same point in space that is *not* on a shared surface.
 *
 * Intersection calculation is done with \c intersect .
 *
 * \note This is a fancier version of the `CellContainer` of GG.
 */
class MaskedUnitTracker final : public UnitTracker
{
  public:
    //@{
    //! Public type aliases
    //@}

  public:
    //// CONSTRUCTION ////

    // Construct from surfaces and regions
    MaskedUnitTracker(SurfaceContainer               surfaces,
                      const std::vector<UnitRegion>& regions);

    //// TRACKING ////

    // Find the local cell and possibly surface ID.
    Initialization initialize(LocalState state) const final;

    // Calculate distance-to-intercept for the next surface
    Intersection intersect(LocalState state) const final;

    // Calculate normal on the current surface
    Real3 normal(LocalState state) const final;

    //// ACCESSORS ////

    //! Number of cells
    size_type num_volumes() const final { return storage_.num_volumes(); }

    //! Number of surfaces
    size_type num_surfaces() const final { return storage_.num_surfaces(); }

    // Access all surfaces
    const SurfaceContainer& surfaces() const final;

    // Reconstruct a CSG cell for diagnostic output
    UnitRegion get_region(VolumeId id) const final;

    // Get the face ID of a surface in a cell
    FaceId find_face(SurfaceId, VolumeId) const final;

  private:
    //// DATA ////

    // Surfaces and cells
    detail::UnitTrackingStorage storage_;
    // Groups of surfaces sorted by connecting zorder
    detail::SurfaceGroupContainer surface_groups_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
