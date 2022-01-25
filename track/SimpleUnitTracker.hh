//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/SimpleUnitTracker.hh
 * \brief SimpleUnitTracker class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "UnitTracker.hh"

#include "Nemesis/containers/Span.hh"
#include "detail/UnitTrackingStorage.hh"

namespace celeritas
{
class SurfaceContainer;
struct UnitRegion;
//---------------------------------------------------------------------------//
/*!
 * Track a particle in a universe of well-connected cells.
 *
 * The simple unit tracker is based on a set of non-overlapping cells
 * comprised of surfaces. It is a faster but less "user-friendly" version of
 * the masked unit tracker because it requires all cells to be exactly
 * defined by their connected surfaces. It does *not* check for overlapping
 * cells.
 */
class SimpleUnitTracker final : public UnitTracker
{
  public:
    //// CONSTRUCTION ////

    // Construct from surfaces and regions
    SimpleUnitTracker(SurfaceContainer               surfaces,
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
    //// TYPEDEFS ////

    using Storage_t        = detail::UnitTrackingStorage;
    using SpanConstCellInt = Storage_t::SpanConstCellInt;

    //// DATA ////

    detail::UnitTrackingStorage storage_;

    //// IMPLEMENTATION METHODS ////

    Intersection simple_intersect(LocalState state) const;
    Intersection complex_intersect(LocalState state) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
