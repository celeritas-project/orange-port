//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/ArrayTrackerTest.hh
 * \brief ArrayTrackerTest class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "TrackerTest.hh"

#include "Nemesis/gtest/GtestFunctions.hh"
#include "base/Casts.hh"
#include "base/Face.hh"
#include "base/VectorFunctions.hh"

namespace orange_test
{
//---------------------------------------------------------------------------//
/*!
 * Helper base class for tests of array trackers.
 *
 * This assumes that:
 * - The tracker encodes the cell as a 3-vector of integer types
 * - inside is the sense when initializing a particle, and outside is the "next
 *   intersection" sense
 * - SurfaceId is a re-encoded Tracker::Face (e.g. HexFace)
 */
//---------------------------------------------------------------------------//

template<class AT>
class ArrayTrackerTest : public TrackerTest
{
    using Base = TrackerTest;

  public:
    //@{
    //! Public type aliases
    using Tracker_t = AT;
    using DimVector = typename Tracker_t::DimVector;
    using Face_t    = typename Tracker_t::Face;
    //@}

  public:
    const Tracker_t& get_tracker() const
    {
        CELER_EXPECT(this->tracker);
        return smart_cast<const Tracker_t&>(*this->tracker);
    }

    VolumeId volume_id(DimVector c) const
    {
        CELER_EXPECT(this->tracker);
        const auto& tracker = this->get_tracker();
        return VolumeId(tracker.volume_id(c));
    }

    virtual SurfaceId surface_id(Face_t f) const
    {
        CELER_EXPECT(f);
        return SurfaceId(f.get());
    }

    using Base::set_state;

    //! Initialize exiting the cell with a given pos/dir
    void set_state(const Real3& pos, const Real3& dir, DimVector cell)
    {
        VolumeId volume_id = this->volume_id(cell);
        return this->set_state(pos, dir, volume_id, {}, {});
    }

    //! Initialize, by default "leaving" the cell
    void set_state(const Real3& pos,
                   const Real3& dir,
                   DimVector    cell,
                   Face_t       face,
                   Sense        sense = Sense::outside)
    {
        CELER_EXPECT(face);
        return this->set_state(
            pos, dir, this->volume_id(cell), this->surface_id(face), sense);
    }

    // Convenience struct for intersect testing
    struct ExpectedIntersection
    {
        // Start direction
        Real3 dir;

        // Expected face/distance
        Face_t    face;
        real_type distance;

        // Expected next cell
        DimVector cell;
    };

    // Convenience struct for intersect testing
    // Test face, distance, and cell for particle initialized at pos
    // in cell with given permutation of direction
    void test_intersect(const Real3&                     pos,
                        const DimVector&                 cell,
                        span<const ExpectedIntersection> tests);
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
template<class AT>
void ArrayTrackerTest<AT>::test_intersect(const Real3&     pos,
                                          const DimVector& cell,
                                          span<const ExpectedIntersection> tests)
{
    const auto& state_ref = this->state_ref();
    for (const auto& t : tests)
    {
        SCOPED_TRACE(make_fixed_view(t.dir));
        this->set_state(pos, t.dir, cell);

        auto init = tracker->initialize(state_ref);
        ASSERT_TRUE(init);
        EXPECT_EQ(Sense::inside, init.sense);

        auto next = tracker->intersect(state_ref);
        ASSERT_TRUE(next);
        EXPECT_SOFT_EQ(t.distance, next.distance);
        EXPECT_EQ(t.face.get(), next.surface.get());
        EXPECT_TRUE(next.surface);
        EXPECT_EQ(Sense::outside, next.sense);

        // Move across surface into next index
        Real3 pos = make_vector(state_ref.pos);
        Real3 dir = make_vector(state_ref.dir);
        axpy(next.distance, dir, pos);
        this->set_state(pos, dir, state_ref.cell, next.surface, next.sense);
        init = tracker->initialize(this->state_ref());
        ASSERT_TRUE(init);
        EXPECT_EQ(this->volume_id(t.cell), init.cell);
    }
}

//---------------------------------------------------------------------------//
} // namespace orange_test

//---------------------------------------------------------------------------//
