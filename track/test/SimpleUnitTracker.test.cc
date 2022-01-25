//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/tstSimpleUnitTracker.cc
 * \brief Tests for class SimpleUnitTracker
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../SimpleUnitTracker.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/Fuzziness.hh"
#include "orange/construct/UnitBuilder.hh"
#include "orange/construct/CuboidShape.hh"
#include "orange/construct/ConeShape.hh"
#include "orange/construct/IntersectionShape.hh"
#include "orange/construct/PlaneShape.hh"
#include "orange/construct/SlabShape.hh"
#include "orange/construct/SphereShape.hh"
#include "UnitTrackerTest.hh"
#include "../detail/Utils.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Utility tests
TEST(CloserFaceTest, all)
{
    celeritas::detail::CloserFace is_closer;
    // Equal
    EXPECT_FALSE(is_closer({FaceId{}, 0.0}, {FaceId{}, 0.0}));
    EXPECT_FALSE(is_closer({FaceId{}, 1.0}, {FaceId{}, 1.0}));

    // Positive vs nonpositive
    EXPECT_FALSE(is_closer({FaceId{}, -0.0001}, {FaceId{}, 1.0}));
    EXPECT_FALSE(is_closer({FaceId{}, 0.0}, {FaceId{}, 1.0}));
    EXPECT_TRUE(is_closer({FaceId{}, 1.0}, {FaceId{}, 0.0}));

    // Positive vs positive
    EXPECT_TRUE(is_closer({FaceId{}, 1.0}, {FaceId{}, 20.0}));
}

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class SimpleUnitTrackerTest : public orange_test::UnitTrackerTest
{
    using Base = orange_test::UnitTrackerTest;

  protected:
    //// TYPE ALIASES ////
    using Tracker_t = celeritas::SimpleUnitTracker;

  protected:
    const Tracker_t& get_tracker() const
    {
        CELER_EXPECT(tracker);
        return dynamic_cast<const Tracker_t&>(*this->tracker);
    }

    using Base::set_state;

    void set_state(const Real3& pos,
                   const Real3& dir,
                   const char*  vol_name,
                   const char*  surf_name,
                   Sense        sense)
    {
        return Base::set_state(pos,
                               dir,
                               this->find_cell(vol_name),
                               this->find_surface(surf_name),
                               sense);
    }

    //! Initialize in the given cell (not on a surface)
    void set_state(const Real3& pos, const Real3& dir, const char* id_to_label)
    {
        return this->set_state(pos, dir, id_to_label, nullptr, Sense::neg);
    }
};

//---------------------------------------------------------------------------//
/*!
 * Single volume, no surfaces.
 *
 * This degenerate case is seen inside a pin cell that's filled with a single
 * material.
 *
 * \todo If/when we have a "filled" unit where anything *not* in a listed
 * volume is in the "fill" cell (this would be good for TRISO matrix and Geant4
 * implementation) we could prohibit this use case, which forces us to relax
 * some requirements.
 */
class InfiniteTest : public SimpleUnitTrackerTest
{
    void SetUp() override
    {
        // Infinite world (clipped implicitly by parent, presumably)
        UnitBuilder build;
        build.region({}, ZOrder::media, ORANGE_MD_FROM_SOURCE("infinite"));
        EXPECT_TRUE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("infinite unit"));
    }
};

TEST_F(InfiniteTest, accessors)
{
    EXPECT_EQ(1, tracker->num_volumes());
    EXPECT_EQ(0, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    std::string expected = R"rst(
======= ======== ============================================================
Cell    Name     Surface logic
======= ======== ============================================================
0       infinite *
                 (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
======= ======== ============================================================

Cells with reentrant surface tracking: "infinite"

======= ==== ============================================================
Surface Name Description
======= ==== ============================================================
======= ==== ============================================================

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(InfiniteTest, initialize_cell)
{
    // Initialize in center
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);
}

TEST_F(InfiniteTest, intersect)
{
    // Initialize in center inside cell 0
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0}, "infinite");
    auto next = tracker->intersect(this->state_ref());
    EXPECT_FALSE(bool(next));
    EXPECT_FALSE(bool(next.surface));
    EXPECT_EQ(no_intersection(), next.distance);
}

TEST_F(InfiniteTest, track)
{
    // Track from an arbitrary position; next distance is always infinite
    auto result = this->track({10, 0, 0}, {0, 0, 1});

    static const char* const expected_cells[] = {"infinite"};
    EXPECT_VEC_EQ(expected_cells, result.cells);
    static const real_type expected_distances[] = {inf};
    EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
}

//---------------------------------------------------------------------------//
/*!
 * Two volumes, one surface.
 *
 * There is a single sphere that divides the universe into two half-spaces,
 * inside and outside.
 */
class SingleSurfaceTest : public SimpleUnitTrackerTest
{
    void SetUp() override
    {
        // Build sphere with radius 1.0 at center (no transform)
        auto sphere = this->shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("sph"), Transform{}, 1.0);

        // Two regions
        UnitBuilder build;
        build.region(
            {{neg, sphere}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("inside"));
        build.region(
            {{pos, sphere}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("outside"));
        EXPECT_TRUE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("single surface"));
    }
};

TEST_F(SingleSurfaceTest, accessors)
{
    EXPECT_EQ(2, tracker->num_volumes());
    EXPECT_EQ(1, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    std::string expected = R"rst(
======= ======= ============================================================
Cell    Name    Surface logic
======= ======= ============================================================
0       inside  0 ~
                (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
1       outside 0
                (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
======= ======= ============================================================

Cells with reentrant surface tracking: "outside"

======= ===== ============================================================
Surface Name  Description
======= ===== ============================================================
0             Sphere: r=1
        sph.s (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
======= ===== ============================================================

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(SingleSurfaceTest, initialize_init)
{
    // Initialize in center
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize in outer area
    this->set_state({10.0, 0.0, 0.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{1}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize on exact boundary (outward)
    this->set_state({1.0, 0.0, 0.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize on exact boundary (inward)
    this->set_state({1.0, 0.0, 0.0}, {-1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);
}

TEST_F(SingleSurfaceTest, initialize_crossing)
{
    // Initialize leaving center cell, as though we've updated the 'next'
    // surface with our post-crossing sense (pos)
    this->set_state({1.0, 0.0, 0.0}, {1, 0, 0}, "inside", "sph.s", pos);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{1}, init.cell);
    EXPECT_EQ(SurfaceId{0}, init.surface);
    EXPECT_EQ(pos, init.sense);

    // Initialize in outer area, headed in
    this->set_state({1.0, 0.0, 0.0}, {-1, 0, 0}, "outside", "sph.s", neg);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{0}, init.surface);
    EXPECT_EQ(neg, init.sense);
}

TEST_F(SingleSurfaceTest, intersect)
{
    // Start in center inside interior
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0}, "inside");
    auto next = tracker->intersect(this->state_ref());
    EXPECT_TRUE(bool(next));
    EXPECT_EQ(SurfaceId{0}, next.surface);
    EXPECT_EQ(pos, next.sense);
    EXPECT_SOFT_EQ(1.0, next.distance);

    // Start on boundary heading out
    this->set_state({1.0, 0.0, 0.0}, {1, 0, 0}, "outside", "sph.s", pos);
    next = tracker->intersect(this->state_ref());
    EXPECT_FALSE(bool(next));
    EXPECT_FALSE(bool(next.surface));
    EXPECT_EQ(no_intersection(), next.distance);

    // Start on outside surface of sphere, heading in
    this->set_state({-1.0, 0.0, 0.0}, {1, 0, 0}, "outside", "sph.s", pos);
    next = tracker->intersect(this->state_ref());
    EXPECT_TRUE(bool(next));
    EXPECT_EQ(SurfaceId{0}, next.surface);
    EXPECT_EQ(neg, next.sense);
#if 0
    EXPECT_SOFT_EQ(1e-8, next.distance); // bump distance
#else
    // Not correct: this skips the zero distance-to-surface because our
    // distance-to-surface routines only account for "on" or "off", not sense.
    // However, this should be mitigated by higher-level logic in the
    // TrackingGeometry that entirely prevents initialization on a lower-level
    // surface.
    EXPECT_SOFT_EQ(2.0, next.distance);
#endif

    // Start on inside surface of sphere, heading in
    this->set_state({-1.0, 0.0, 0.0}, {1, 0, 0}, "inside", "sph.s", neg);
    next = tracker->intersect(this->state_ref());
    EXPECT_TRUE(bool(next));
    EXPECT_EQ(SurfaceId{0}, next.surface);
    EXPECT_EQ(pos, next.sense);
    EXPECT_SOFT_EQ(2.0, next.distance); // bump distance

    // Start exactly on the surface, but without a surface/sense. This
    // can happen when we're tracking in a lower universe but are calculating
    // the next distance at a higher level.
    this->set_state({-1.0, 0.0, 0.0}, {1, 0, 0}, "inside");
    next = tracker->intersect(this->state_ref());
    EXPECT_TRUE(bool(next));
    EXPECT_EQ(SurfaceId{0}, next.surface);
    EXPECT_EQ(pos, next.sense);
    EXPECT_SOFT_EQ(2.0, next.distance);

    this->set_state({1.0, 0.0, 0.0}, {1, 0, 0}, "inside");
    next = tracker->intersect(this->state_ref());
    // TODO: this behavior is currently incorrect: we should be able to give an
    // exiting distance of zero (epsilon) without registering as being on the
    // surface. We might just have to have "on surface/sense" for every level
    // to really do this right, and that might also allows us to do
    // zero-distance
#if 0
    EXPECT_TRUE(bool(next));
    EXPECT_EQ(SurfaceId{0}, next.surface);
    EXPECT_EQ(pos, next.sense);
    EXPECT_EQ(0.0, next.distance);
#else
    EXPECT_FALSE(bool(next));
    EXPECT_FALSE(bool(next.surface));
    EXPECT_EQ(no_intersection(), next.distance);
#endif
}

TEST_F(SingleSurfaceTest, track_leftward)
{
    auto result = this->track({10, 0, 0}, {-1, 0, 0});

    static const char* const expected_cells[]
        = {"outside", "inside", "outside"};
    EXPECT_VEC_EQ(expected_cells, result.cells);
    static const char* const expected_surfaces[] = {"", "sph.s", "sph.s"};
    EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
    static const int expected_senses[] = {' ', '-', '+'};
    EXPECT_VEC_EQ(expected_senses, result.senses);
    static const real_type expected_distances[] = {9, 2, inf};
    EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
}

TEST_F(SingleSurfaceTest, track_from_initial_coincident)
{
    auto result = this->track({-1, 0, 0}, {1, 0, 0});

    static const char* const expected_cells[] = {"[fail]"};
    EXPECT_VEC_EQ(expected_cells, result.cells);
}

//---------------------------------------------------------------------------//
/*!
 * Plane edge initialization.
 *
 * There is an 'inside' and outside.
 */
class PlaneEdgeTest : public SimpleUnitTrackerTest
{
    void SetUp() override
    {
        using PairDbl = ConeShape::PairDbl;

        // clang-format off
        auto slab = this->shapes.emplace<SlabShape>(
            ORANGE_MD_FROM_SOURCE("global"), Transform{},
            Axis::x, -1.0, 1.0);
        // clang-format on

        UnitBuilder build;
        build.region(
            {{pos, slab}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("exterior"));
        build.region(
            {{neg, slab}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("inside"));
        EXPECT_TRUE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("planeedge"));
    }
};

// Test initialization edge cases
TEST_F(PlaneEdgeTest, track)
{
    {
        SCOPED_TRACE("rightward starting outside");
        auto result = this->track({-1.5, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "inside", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {0.5, 2, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        SCOPED_TRACE("rightward coincident");
        auto result = this->track({-1., 0, 0}, {1, 0, 0});

        static const char* const expected_cells[] = {"[fail]"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
    }
    {
        SCOPED_TRACE("leftward coincident");
        auto                     result = this->track({1., 0, 0}, {-1, 0, 0});
        static const char* const expected_cells[] = {"[fail]"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Interior cone.
 *
 * The left cell (a cuboid) has an interior cone with one surface coincident
 * with the right boundary. There is also a cell on the right side opposite the
 * cone.
 */
class InteriorConeTest : public SimpleUnitTrackerTest
{
    void SetUp() override
    {
        using PairDbl = ConeShape::PairDbl;

        // clang-format off
        auto left = this->shapes.emplace<CuboidShape>(
            ORANGE_MD_FROM_SOURCE("left"), Transform{},
            Real3{0, 0, -1}, Real3{8, 6, 1});
        auto cone = this->shapes.emplace<ConeShape>(
            ORANGE_MD_FROM_SOURCE("cone"), Transform{Real3{0, 3, 0}},
            Axis::x, PairDbl{0, 2}, PairDbl{4, 8});
        auto right = this->shapes.emplace<CuboidShape>(
            ORANGE_MD_FROM_SOURCE("right"), Transform{},
            Real3{8, 0, -1}, Real3{10, 6, 1});
        auto global = this->shapes.emplace<CuboidShape>(
            ORANGE_MD_FROM_SOURCE("global"), Transform{},
            Real3{0, 0, -1}, Real3{10, 6, 1});
        // clang-format on

        UnitBuilder build;
        build.region(
            {{pos, global}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("exterior"));
        build.region(
            {{neg, cone}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("cone"));
        build.region({{neg, left}, {pos, cone}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("left"));
        build.region(
            {{neg, right}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("right"));
        EXPECT_TRUE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("interiorcone"));
    }
};

TEST_F(InteriorConeTest, accessors)
{
    EXPECT_EQ(4, tracker->num_volumes());
    EXPECT_EQ(6 + 3, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    std::string expected = R"rst(
======= ======== ============================================================
Cell    Name     Surface logic
======= ======== ============================================================
0       exterior 0 1 ~ & 2 & 3 ~ & 4 & 5 ~ & ~
                 (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
1       cone     6 ~ 7 & 8 ~ &
                 (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
2       left     6 ~ 7 & 8 ~ & ~ 0 2 & 3 ~ & 4 & 5 ~ & 8 ~ & &
                 (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
3       right    1 ~ 2 & 3 ~ & 4 & 5 ~ & 8 &
                 (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
======= ======== ============================================================

Cells with reentrant surface tracking: "exterior", "left"

======= ========= ============================================================
Surface Name      Description
======= ========= ============================================================
0                 Plane: x=0
        global.mx (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        left.mx   (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
1                 Plane: x=10
        global.px (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        right.px  (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
2                 Plane: y=0
        global.my (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        left.my   (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        right.my  (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
3                 Plane: y=6
        global.py (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        left.py   (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        right.py  (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
4                 Plane: z=-1
        global.mz (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        left.mz   (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        right.mz  (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
5                 Plane: z=1
        global.pz (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        left.pz   (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        right.pz  (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
6                 Cone x: tangent=0.5 at 4 3 0
        cone.kox  (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
7                 Plane: x=4
        cone.mx   (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
8                 Plane: x=8
        cone.px   (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        left.px   (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
        right.mx  (from ``track/test/tstSimpleUnitTracker.cc:NNN``)
======= ========= ============================================================

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(InteriorConeTest, move_leftward)
{
    // Initialize leaving exterior cell into right cell
    this->set_state({11.0, 0.5, 0.0}, {-1, 0, 0});
    EXPECT_EQ(VolumeId{}, this->state_ref().cell);
    EXPECT_EQ(SurfaceId{}, this->state_ref().surface);

    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize leaving exterior cell into right cell
    this->set_state({10.0, 0.5, 0.0}, {-1, 0, 0}, "exterior", "global.px", neg);
    EXPECT_EQ(VolumeId{0}, this->state_ref().cell);
    EXPECT_EQ(SurfaceId{1}, this->state_ref().surface);
    EXPECT_EQ(Sense::neg, this->state_ref().sense);

    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{3}, init.cell);
    EXPECT_EQ(SurfaceId{1}, init.surface);
    EXPECT_EQ(Sense::neg, this->state_ref().sense);
}

TEST_F(InteriorConeTest, track)
{
    {
        // Track rightward, not hitting cone
        auto result = this->track({-1, 0.5, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "left", "right", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "cone.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 8, 2, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track rightward, through cone tangent point
        auto result = this->track({-1, 4, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "left", "cone", "right", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "cone.kox", "cone.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '-', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 6, 2, 2, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track leftward, missing cone
        auto result = this->track({11, 0.5, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "right", "left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.px", "cone.px", "global.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 2, 8, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track up across cone surface but without entering the cone cell
        auto result = this->track({1, -7, 0}, {0, 1, 0});

        static const char* const expected_cells[]
            = {"exterior", "left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.my", "global.py"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {7, 6, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track up, actually going through cone
        auto result = this->track({7, -7, 0}, {0, 1, 0});

        static const char* const expected_cells[]
            = {"exterior", "left", "cone", "left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.my", "cone.kox", "cone.kox", "global.py"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '-', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {7, 1.5, 3, 1.5, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}
