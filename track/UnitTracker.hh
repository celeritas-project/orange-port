//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/UnitTracker.hh
 * \brief UnitTracker class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Tracker.hh"

namespace celeritas
{
class SurfaceContainer;
struct UnitRegion;
//---------------------------------------------------------------------------//
/*!
 * Interface class for CSG trackers
 *
 * This facilitates 'metadata' output for any low-level unit tracker
 * implementation.
 */
class UnitTracker : public Tracker
{
  public:
    //! Access all surfaces
    virtual const SurfaceContainer& surfaces() const = 0;

    //! Reconstruct a CSG cell for diagnostic output
    virtual UnitRegion get_region(VolumeId id) const = 0;

    //! Get the face ID of a surface in a cell
    virtual FaceId find_face(SurfaceId, VolumeId) const = 0;

    //! Exterior cell is always the zeroth
    static constexpr VolumeId exterior_volume() { return VolumeId{0}; }

  protected:
    // Don't allow deletion via UnitTracker pointer: we're just an interface!
    ~UnitTracker() = default;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
