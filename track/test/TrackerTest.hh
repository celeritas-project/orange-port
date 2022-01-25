//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/TrackerTest.hh
 * \brief TrackerTest class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Nemesis/gtest/Test.hh"

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "orange/Definitions.hh"
#include "../Definitions.hh"
#include "../SenseContainer.hh"

namespace celeritas
{
class Tracker;
class UniverseMetadata;
} // namespace celeritas

namespace orange_test
{
//---------------------------------------------------------------------------//
/*!
 * Helper class for testing tracker functionality.
 */
class TrackerTest : public ::Test
{
  public:
    //@{
    //! Public type aliases
    using Real3          = celeritas::Real3;
    using VolumeId       = celeritas::VolumeId;
    using SurfaceId      = celeritas::SurfaceId;
    using Sense          = celeritas::Sense;
    using Initialization = celeritas::Initialization;
    using Intersection   = celeritas::Intersection;

    using LocalState = celeritas::LocalState;
    using UPConstTracker     = std::unique_ptr<const celeritas::Tracker>;
    using SPConstUnivMetadata
        = std::shared_ptr<const celeritas::UniverseMetadata>;
    //@}

    struct TrackResult
    {
        std::vector<std::string> cells;
        std::vector<std::string> surfaces;
        std::vector<char>        senses;
        std::vector<real_type>   distances;

        void print_expected() const;
    };

#ifndef ORANGE_INTEL_CONSTEXPR_BUG
    static constexpr real_type inf = celeritas::no_intersection();
#else
    static const real_type inf;
#endif

  public:
    // Constructor
    TrackerTest();

    //// TRACKING ////

    // Initialize with a given pos/dir
    void set_state(const Real3& pos,
                   const Real3& dir,
                   VolumeId     cell,
                   SurfaceId    surface,
                   Sense        sense);

    // Initialize with a given pos/dir
    void set_state(const Real3& pos, const Real3& dir)
    {
        return this->set_state(pos, dir, {}, {}, Sense::inside);
    }

    //! Get a state reference
    const LocalState& state_ref() const { return state_ref_; }

    // Track until an infinite distance/exception/fixed point is reached
    TrackResult track(const Real3& pos, const Real3& dir);

    //// QUERYING ////

    // Find the cell from its label (nullptr allowed)
    VolumeId find_cell(const char* label) const;

    // Find the surface from its label (NULL pointer allowed)
    SurfaceId find_surface(const char* label) const;

    // Surface name (or sentinel if no surface);
    std::string id_to_label(SurfaceId surf) const;

    // Cell name (or sentinel if no surface);
    std::string id_to_label(VolumeId) const;

    // Get the string output from the metadata
    std::string describe_md() const;

  protected:
    void set_tracker(UPConstTracker tracker, SPConstUnivMetadata md);

    //// DATA (read-only please!) ////

    UPConstTracker      tracker;
    SPConstUnivMetadata md;

  private:
    Real3                     pos_;
    Real3                     dir_;
    celeritas::SenseContainer temp_senses_;
    celeritas::VecNextFace     temp_face_dist_;
    LocalState        state_ref_;

    std::map<std::string, VolumeId>  volume_ids_;
    std::map<std::string, SurfaceId> surface_ids_;
};

//---------------------------------------------------------------------------//
} // namespace orange_test

//---------------------------------------------------------------------------//
