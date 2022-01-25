//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/Universe.hh
 * \brief Universe class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "orange/Transform.hh"
#include "TrackingState.hh"
#include "Tracker.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Dispatch to Trackers and associate cells with daughters.
 *
 * The Universe is a very thin class that manages daughter universes and
 * dispatches to the stored Tracker class. It is the *only* tracking class that
 * should contain a virtual function call, and can thus be used to hide
 * the tracker's polymorphism implementation if needed.
 *
 * A universe will dispatch these to its associated tracker:
 * - Initialization
 * - Surface crossing
 * - Distance-to-boundary calculation
 *
 * The universe does *not* do any movement operations, nor does it modify the
 * state (except for its use of temporary state storage).
 *
 * The daughter array implementation is optimized for fast construction and
 * access, targeting the worst cases of "all local cells and no universes" (a
 * complicated lowest-level unit) and "no local cells and all universes" (a
 * very large array).
 */
class Universe
{
  public:
    //@{
    //! Public type aliases
    using State_t         = TrackingState;
    using size_type       = VolumeId::size_type;
    using SPConstUniverse = std::shared_ptr<const Universe>;
    using UPConstTracker  = std::unique_ptr<const Tracker>;
    //@}

    struct Daughter
    {
        SPConstUniverse universe;
        Transform       transform; //!< Daughter-to-parent
    };

    using VecDaughter = std::vector<Daughter>;

    //! Construction parameters
    struct Params
    {
        UniverseId     id;
        UPConstTracker tracker;
        VecDaughter    daughters; //!< Empty or one per cell
    };

  public:
    // Construct with parameters
    explicit Universe(Params params);

    // Initialize a particle or cross a surface
    inline Initialization initialize(LocalState) const;

    // Find the distance-to-intercept at this level
    inline Intersection intersect(LocalState) const;

    // Calculate the *local* surface normal
    inline Real3 normal(LocalState) const;

    //// ACCESSORS ////

    //! Universe ID
    UniverseId id() const { return id_; }

    //! Access the tracker (primarily for output purposes?)
    const Tracker& tracker() const { return *tracker_; }

    //! Number of cells in this universe
    size_type num_volumes() const { return tracker_->num_volumes(); }

    //! Number of surfaces in this universe
    size_type num_surfaces() const { return tracker_->num_surfaces(); }

    // Get the local daughter we're currently in; nullptr if none
    inline const Daughter* daughter(VolumeId cell) const;

  private:
    UniverseId     id_;
    UPConstTracker tracker_;
    VecDaughter    daughters_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Universe.i.hh"
//---------------------------------------------------------------------------//
