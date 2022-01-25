//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/TrackingState.hh
 * \brief TrackingState class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Macros.hh"
#include "Definitions.hh"
#include "SenseContainer.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * \struct TrackingState
 * State of a particle use for tracking.
 *
 * The tracking state only cares about the geometry, not materials or labels or
 * anything like that.
 *
 * TODO: use absl::InlinedVector or equivalent instead of fixed-sized max
 * level. As a point of comparison, VecGeom precalculates "max level" from the
 * geometry and can allocate a bunch of states simultaneously in a "state
 * pool". When porting to GPU  it might be feasible to have a
 * `StateAllocation`
 * class that contains a device vector of `LevelState` and another device
 * vector of `TrackingState`, and change `level_state` to a pointer that
 * references `StateAllocation`-owned data.
 */
//---------------------------------------------------------------------------//

struct TrackingState
{
    //@{
    //! Type aliases
    using NextFace      = std::pair<FaceId, real_type>;
    using VecNextFace   = std::vector<NextFace>;
    using Intersection = ::celeritas::Intersection;
    //@}

    //! State local to a level
    struct LevelState
    {
        Real3        pos;
        Real3        dir;
        VolumeId     cell;
        UniverseId   univ_id;
        Intersection next;

        template<class Archiver>
        void serialize(Archiver& ar)
        {
            // clang-format off
            ar & pos & dir & cell & univ_id & next;
            // clang-format on
        }
    };

    struct OnSurface
    {
        int       level = -1;
        SurfaceId surface;
        Sense     sense;

        //! Whether we're on a surface
        explicit operator bool() const { return bool(surface); }

        template<class archiver>
        void serialize(archiver& ar)
        {
            // clang-format off
            ar & level & surface & sense;
            // clang-format on
        }
    };

    //! Maximum number of levels deep
    static constexpr int max_level = 20;

    //! Current level
    int level = -1;

    //! Level of next intersection
    int next_level = -1;

    //! State at every level
    LevelState level_state[max_level];

    //! Surface IDs in the current level if on a surface
    OnSurface on_surface;

    //! Cumulative straight-line distance traveled since last surface crossing
    //! or direction change
    real_type movement = 0;

    //! Cumulative error term from bumping (add to pos to get "actual")
    Real3 correction;

    //@{
    //! Scratch space
    mutable SenseContainer temp_senses;
    mutable VecNextFace     temp_face_dist;
    //@}

    //! Whether the state is initialized
    CELER_FORCEINLINE_FUNCTION explicit operator bool() const
    {
        return this->level >= 0;
    }

    template<class Archiver>
    void serialize(Archiver& ar)
    {
        // clang-format off
        ar & level & next_level & on_surface & movement & correction;
        // clang-format on

        int i = 0;
        for (; i <= level; ++i)
            ar& level_state[i];

        // Since packed states must have fixed size, zip up a dummy level for
        // remaining levels.
        LevelState dummy;
        for (; i < max_level; ++i)
            ar& dummy;
    }
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
