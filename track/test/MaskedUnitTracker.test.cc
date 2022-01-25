//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/tstMaskedUnitTracker.cc
 * \brief Tests for class MaskedUnitTracker
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#include "../MaskedUnitTracker.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/Fuzziness.hh"
#include "orange/construct/UnitBuilder.hh"
#include "orange/construct/CuboidShape.hh"
#include "orange/construct/ConeShape.hh"
#include "orange/construct/IntersectionShape.hh"
#include "orange/construct/PlaneShape.hh"
#include "orange/construct/SphereShape.hh"
#include "UnitTrackerTest.hh"

using namespace celeritas;
using constants::sqrt_two;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class MaskedUnitTrackerTest : public orange_test::UnitTrackerTest
{
    using Base = orange_test::UnitTrackerTest;

  protected:
    //// TYPE ALIASES ////
    using Tracker_t = celeritas::MaskedUnitTracker;

  protected:
    const Tracker_t& get_tracker() const
    {
        CELER_EXPECT(tracker);
        return dynamic_cast<const Tracker_t&>(*this->tracker);
    }

    // Initialize by finding a cell
    using Base::set_state;

    // Initialize in a cell
    void set_state(const Real3& pos, const Real3& dir, const char* id_to_label)
    {
        return Base::set_state(
            pos, dir, this->find_cell(id_to_label), {}, Sense::neg);
    }

    // Initialize on a surface
    void set_state(const Real3& pos,
                   const Real3& dir,
                   const char*  id_to_label,
                   const char*  id_to_label,
                   Sense        sense)
    {
        return Base::set_state(pos,
                               dir,
                               this->find_cell(id_to_label),
                               this->find_surface(id_to_label),
                               sense);
    }
};

//---------------------------------------------------------------------------//
/*!
 * Single volume, no surfaces.
 *
 * In practice we should probably forbid this (since it could be replaced by
 * the simple tracker) and can remove this test. However, it's useful for unit
 * testing.
 */
class InfiniteTest : public MaskedUnitTrackerTest
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

//---------------------------------------------------------------------------//
/*!
 * Two volumes, one surface.
 *
 * This also can be replaced by a simple unit but should nevertheless be
 * correct here.
 */
class SingleSurfaceTest : public MaskedUnitTrackerTest
{
    void SetUp() override
    {
        // Build sphere with radius 1.0 at center (no transform)
        auto sphere = this->shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("sph"), Transform{}, 1.0);

        // Two regions
        UnitBuilder build;
        build.region(
            {{inside, sphere}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("inside"));
        build.region({{outside, sphere}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("outside"));
        EXPECT_TRUE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("single surface"));
    }
};

TEST_F(SingleSurfaceTest, accessors)
{
    EXPECT_EQ(2, tracker->num_volumes());
    EXPECT_EQ(1, tracker->num_surfaces());
}

TEST_F(SingleSurfaceTest, initialize_cell)
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
}

TEST_F(SingleSurfaceTest, initialize_surface)
{
    // Initialize leaving center cell
    this->set_state({1.0, 0.0, 0.0}, {1, 0, 0}, "inside", "sph.s", pos);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{1}, init.cell);
    EXPECT_EQ(SurfaceId{0}, init.surface);
    EXPECT_EQ(Sense::pos, init.sense);

    // Initialize in outer area, headed in
    this->set_state({1.0, 0.0, 0.0}, {-1, 0, 0}, "outside", "sph.s", neg);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{0}, init.surface);
    EXPECT_EQ(Sense::neg, init.sense);
}

TEST_F(SingleSurfaceTest, intersect)
{
    // Start in center inside interior
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0}, "inside");
    auto next = tracker->intersect(this->state_ref());
    EXPECT_TRUE(bool(next));
    EXPECT_EQ(SurfaceId{0}, next.surface);
    EXPECT_EQ(Sense::pos, next.sense);
    EXPECT_EQ(1.0, next.distance);

    // Start on boundary heading out
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0}, "outside", "sph.s", pos);
    next = tracker->intersect(this->state_ref());
    EXPECT_FALSE(bool(next));
    EXPECT_FALSE(bool(next.sense));
    EXPECT_FALSE(bool(next.surface));
    EXPECT_EQ(no_intersection(), next.distance);
}

TEST_F(SingleSurfaceTest, track)
{
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0}, "inside");
}

//---------------------------------------------------------------------------//
/*!
 * Two volumes, one is a hole.
 *
 * A single sphere "hole" masks out the infinite background medium.
 */
class SingleHoleTest : public MaskedUnitTrackerTest
{
  public:
    using PlaneVector = Array<real_type, 2>;

    void SetUp() override
    {
        // Build sphere with radius 1.0 at center (no transform)
        auto sphere = this->shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("sph"), Transform{}, 1.0);

        // Two regions
        UnitBuilder build;
        build.region(
            {{inside, sphere}}, ZOrder::hole, ORANGE_MD_FROM_SOURCE("hole"));
        build.region({}, ZOrder::media, ORANGE_MD_FROM_SOURCE("everywhere"));
        EXPECT_FALSE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("single hole"));
    }
};

TEST_F(SingleHoleTest, accessors)
{
    EXPECT_EQ(2, tracker->num_volumes());
    EXPECT_EQ(1, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    std::string expected = R"rst(
======= = ========== ============================================================
Cell    Z Name       Surface logic
======= = ========== ============================================================
0       H hole       0 ~
                     (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
1       M everywhere *
                     (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
======= = ========== ============================================================

Cells with reentrant surface tracking: "everywhere"

======= ===== ============================================================
Surface Name  Description
======= ===== ============================================================
0             Sphere: r=1
        sph.s (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
======= ===== ============================================================

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(SingleHoleTest, initialize_cell)
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
}

TEST_F(SingleHoleTest, initialize_surface)
{
    // Initialize leaving center cell
    this->set_state({1.0, 0.0, 0.0}, {1, 0, 0}, "hole", "sph.s", pos);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{1}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize in outer area, headed in
    this->set_state({1.0, 0.0, 0.0}, {-1, 0, 0}, "everywhere", "sph.s", neg);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{0}, init.surface);
    EXPECT_EQ(Sense::neg, init.sense);
}

TEST_F(SingleHoleTest, intersect)
{
    // Start in center inside interior
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0}, "hole");
    auto next = tracker->intersect(this->state_ref());
    EXPECT_TRUE(bool(next));
    EXPECT_EQ(SurfaceId{0}, next.surface);
    EXPECT_EQ(Sense::pos, next.sense);
    EXPECT_EQ(1.0, next.distance);

    // Start on boundary heading out
    this->set_state({0.0, 0.0, 0.0}, {1, 0, 0}, "everywhere", "sph.s", pos);
    next = tracker->intersect(this->state_ref());
    EXPECT_FALSE(bool(next));
    EXPECT_FALSE(bool(next.sense));
    EXPECT_FALSE(bool(next.surface));
    EXPECT_EQ(no_intersection(), next.distance);
}

//---------------------------------------------------------------------------//
/*!
 * INTERIOR CONE (cone is a 'hole' here but tracking should be identical to
 * simple unit, except for coincident faces when entering/leaving hole).
 */
class InteriorConeTest : public MaskedUnitTrackerTest
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
        build.region({{outside, global}},
                     ZOrder::exterior,
                     ORANGE_MD_FROM_SOURCE("exterior"));
        build.region(
            {{inside, cone}}, ZOrder::hole, ORANGE_MD_FROM_SOURCE("cone"));
        build.region(
            {{inside, left}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("left"));
        build.region(
            {{inside, right}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("right"));
        EXPECT_FALSE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("interiorcone"));
    }
};

TEST_F(InteriorConeTest, accessors)
{
    EXPECT_EQ(4, tracker->num_volumes());
    EXPECT_EQ(9, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);
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
    EXPECT_EQ(Sense::neg, init.sense);
}

TEST_F(InteriorConeTest, track_horizontal)
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
        // Track reversed, missing cone
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
        // Track rightward, through cone tangent point
        auto result = this->track({-1, 3, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "left", "cone", "right", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "cone.kox", "cone.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '-', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 4, 4, 2, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track reversed (starting 1 unit outside again)
        auto result = this->track({11, 3, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "right", "cone", "left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.px", "cone.px", "", "global.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-', '-', ' ', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 2, 4, 4, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

TEST_F(InteriorConeTest, track_vertical)
{
    {
        // Track up across cone surface but without entering the cone cell
        auto result = this->track({1, -1, 0}, {0, 1, 0});

        static const char* const expected_cells[]
            = {"exterior", "left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.my", "global.py"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 6, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track up, actually going through cone
        auto result = this->track({7, -1, 0}, {0, 1, 0});

        static const char* const expected_cells[]
            = {"exterior", "left", "cone", "left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.my", "cone.kox", "", "global.py"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '-', ' ', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1.5, 3, 1.5, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Square hole.
 */
class SquareHoleTest : public MaskedUnitTrackerTest
{
    void SetUp() override
    {
        auto outer
            = this->shapes.emplace<CuboidShape>(ORANGE_MD_FROM_SOURCE("outer"),
                                                Transform{},
                                                Real3{0, 0, -1},
                                                Real3{10, 10, 1});
        auto up_left = this->shapes.emplace<CuboidShape>(
            ORANGE_MD_FROM_SOURCE("upper left"),
            Transform{},
            Real3{1, 5, -1},
            Real3{5, 9, 1});
        auto lo_right = this->shapes.emplace<CuboidShape>(
            ORANGE_MD_FROM_SOURCE("lower left"),
            Transform{},
            Real3{5, 1, -1},
            Real3{9, 5, 1});
        auto center = this->shapes.emplace<CuboidShape>(
            ORANGE_MD_FROM_SOURCE("center"),
            Transform{},
            Real3{3, 3, -1},
            Real3{7, 7, 1});

        // Two regions
        UnitBuilder build;
        build.region({{outside, outer}},
                     ZOrder::exterior,
                     ORANGE_MD_FROM_SOURCE("outside"));
        build.region(
            {{inside, center}}, ZOrder::hole, ORANGE_MD_FROM_SOURCE("hole"));
        build.region({{inside, up_left}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("media 1"));
        build.region({{inside, lo_right}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("media 2"));
        build.region({{inside, outer}, {outside, up_left}, {outside, lo_right}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("media 3"));
        EXPECT_FALSE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("square hole"));
    }
};

TEST_F(SquareHoleTest, accessors)
{
    EXPECT_EQ(5, tracker->num_volumes());
    EXPECT_EQ(16, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    std::string expected = R"rst(
======= = ======= ============================================================
Cell    Z Name    Surface logic
======= = ======= ============================================================
0       B outside 0 1 ~ & 2 & 3 ~ & 4 & 5 ~ & ~
                  (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
1       H hole    4 5 ~ & 6 & 7 ~ & 8 & 9 ~ &
                  (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
2       M media 1 4 5 ~ & 10 & 11 ~ & 12 & 13 ~ &
                  (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
3       M media 2 4 5 ~ & 11 & 12 ~ & 14 ~ & 15 &
                  (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
4       M media 3 0 1 ~ & 2 & 3 ~ & 4 & 5 ~ & 4 5 ~ & 10 & 11 ~ & 12 & 13 ~ &
                  ~ & 4 5 ~ & 11 & 12 ~ & 14 ~ & 15 & ~ &
                  (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
======= = ======= ============================================================

Cells with reentrant surface tracking: "outside", "media 3"

======= ============= ============================================================
Surface Name          Description
======= ============= ============================================================
0                     Plane: x=0
        outer.mx      (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
1                     Plane: x=10
        outer.px      (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
2                     Plane: y=0
        outer.my      (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
3                     Plane: y=10
        outer.py      (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
4                     Plane: z=-1
        center.mz     (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
        lower left.mz (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
        outer.mz      (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
        upper left.mz (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
5                     Plane: z=1
        center.pz     (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
        lower left.pz (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
        outer.pz      (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
        upper left.pz (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
6                     Plane: x=3
        center.mx     (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
7                     Plane: x=7
        center.px     (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
8                     Plane: y=3
        center.my     (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
9                     Plane: y=7
        center.py     (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
10                    Plane: x=1
        upper left.mx (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
11                    Plane: x=5
        lower left.mx (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
        upper left.px (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
12                    Plane: y=5
        lower left.py (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
        upper left.my (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
13                    Plane: y=9
        upper left.py (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
14                    Plane: x=9
        lower left.px (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
15                    Plane: y=1
        lower left.my (from ``track/test/tstMaskedUnitTracker.cc:NNN``)
======= ============= ============================================================

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(SquareHoleTest, track_horizontal)
{
    {
        // Track rightward, entering hole at (3, 4) and leaving at (7, 4)
        auto result = this->track({-1, 4, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"outside", "media 3", "hole", "media 2", "media 3", "outside"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "outer.mx", "center.mx", "", "lower left.px", "outer.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+', ' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 3, 4, 2, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        auto result = this->track({11, 4, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"outside", "media 3", "media 2", "hole", "media 3", "outside"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "outer.px", "lower left.px", "center.px", "", "outer.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-', '-', '-', ' ', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 2, 4, 3, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

TEST_F(SquareHoleTest, track_diagonal)
{
    const real_type eps = celeritas::fuzziness().bump_abs() * 0.1;

    {
        // Track down and right
        auto result = this->track({0 + eps, 11 - eps, 0}, {1, -1, 0});

        static const char* const expected_cells[] = {"outside",
                                                     "media 3",
                                                     "media 1",
                                                     "hole",
                                                     "media 2",
                                                     "media 3",
                                                     "outside"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {"",
                                                        "outer.py",
                                                        "upper left.py",
                                                        "center.py",
                                                        "",
                                                        "lower left.px",
                                                        "outer.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[]
            = {' ', '-', '-', '-', ' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1.41421356095888,
                                                       1.414213562373,
                                                       2.828427124746,
                                                       4.242640687119,
                                                       2.828427124746,
                                                       1.414213562373,
                                                       inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        auto result = this->track({11 - eps, 0 + eps, 0}, {-1, 1, 0});

        static const char* const expected_cells[] = {"outside",
                                                     "media 3",
                                                     "media 2",
                                                     "hole",
                                                     "media 1",
                                                     "media 3",
                                                     "outside"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {"",
                                                        "outer.px",
                                                        "lower left.px",
                                                        "center.px",
                                                        "",
                                                        "upper left.py",
                                                        "outer.py"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[]
            = {' ', '-', '-', '-', ' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1.41421356095888,
                                                       1.414213562373,
                                                       2.828427124746,
                                                       4.242640687119,
                                                       2.828427124746,
                                                       1.414213562373,
                                                       inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

TEST_F(SquareHoleTest, track_starting_in_hole)
{
    // Track out through hole
    auto result = this->track({4, 6, 0}, {0, 1, 0});

    static const char* const expected_cells[]
        = {"hole", "media 1", "media 3", "outside"};
    EXPECT_VEC_EQ(expected_cells, result.cells);
    static const char* const expected_surfaces[]
        = {"", "", "upper left.py", "outer.py"};
    EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
    static const int expected_senses[] = {' ', ' ', '+', '+'};
    EXPECT_VEC_EQ(expected_senses, result.senses);
    static const real_type expected_distances[] = {1, 2, 1, inf};
    EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
}

TEST_F(SquareHoleTest, cut_corners_tight)
{
    // Skips across corners
    const real_type eps = celeritas::fuzziness().bump_abs() * 0.1;

    {
        // Track down and right, barely inside corners
        auto result = this->track({-1, 7 + eps, 0}, {1, -1, 0});

        static const char* const expected_cells[]
            = {"outside", "media 3", "outside"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "outer.mx", "outer.my"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const real_type expected_distances[]
            = {sqrt_two, (6 + eps) * sqrt_two, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        auto result = this->track({7 + eps, -1, 0}, {-1, 1, 0});
        static const char* const expected_cells[]
            = {"outside", "media 3", "outside"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "outer.my", "outer.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[]
            = {sqrt_two, (6 + eps) * sqrt_two, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

TEST_F(SquareHoleTest, cut_corners_loose)
{
    const real_type eps = celeritas::fuzziness().bump_abs() * 5;
    {
        // Track down and right, barely inside corners
        auto result = this->track({-1, 7 + eps, 0}, {1, -1, 0});

        static const char* const expected_cells[] = {"outside",
                                                     "media 3",
                                                     "media 1",
                                                     "media 3",
                                                     "hole",
                                                     "media 3",
                                                     "media 2",
                                                     "media 3",
                                                     "outside"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {"",
                                                        "outer.mx",
                                                        "upper left.mx",
                                                        "lower left.py",
                                                        "center.mx",
                                                        "",
                                                        "lower left.mx",
                                                        "lower left.my",
                                                        "outer.my"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[]
            = {' ', '+', '+', '-', '+', ' ', '+', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1.414213562373,
                                                       1.414213562373,
                                                       7.071067768891e-08,
                                                       2.828427054036,
                                                       7.071067768891e-08,
                                                       2.828427054036,
                                                       7.071067768891e-08,
                                                       1.414213562373,
                                                       inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Reversed
        auto result = this->track({7 + eps, -1, 0}, {-1, 1, 0});

        static const char* const expected_cells[] = {"outside",
                                                     "media 3",
                                                     "media 2",
                                                     "media 3",
                                                     "hole",
                                                     "media 3",
                                                     "media 1",
                                                     "media 3",
                                                     "outside"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {"",
                                                        "outer.my",
                                                        "lower left.my",
                                                        "lower left.mx",
                                                        "center.my",
                                                        "",
                                                        "lower left.py",
                                                        "upper left.mx",
                                                        "outer.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[]
            = {' ', '+', '+', '-', '+', ' ', '+', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1.414213562373,
                                                       1.414213562373,
                                                       7.071067768891e-08,
                                                       2.828427054036,
                                                       7.071067768891e-08,
                                                       2.828427054036,
                                                       7.071067768891e-08,
                                                       1.414213562373,
                                                       inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Coincident squares: media is divided into 4 quadrants, and there's an
 * arbitrary shape that overlaps part of the problem.
 */
class CoincidentSquareTest : public MaskedUnitTrackerTest
{
  public:
    using PlaneVector = Array<real_type, 2>;
    using VecInt      = std::vector<int>;

    SPConstShape make_square(ObjectMetadata md, PlaneVector lo, PlaneVector hi)
    {
        return this->shapes.emplace<CuboidShape>(std::move(md),
                                                 Transform{},
                                                 Real3{lo[0], lo[1], -1},
                                                 Real3{hi[0], hi[1], 1});
    }

    void build(SPConstShape hole_shape)
    {
        using PairDbl = ConeShape::PairDbl;

        auto upright = this->make_square(
            ORANGE_MD_FROM_SOURCE("upper-right"), {0, 0}, {1, 1});
        auto upleft = this->make_square(
            ORANGE_MD_FROM_SOURCE("upper-left"), {-1, 0}, {0, 1});
        auto loleft = this->make_square(
            ORANGE_MD_FROM_SOURCE("lower-left"), {-1, -1}, {0, 0});
        auto loright = this->make_square(
            ORANGE_MD_FROM_SOURCE("lower-right"), {0, -1}, {1, 0});

        auto global = this->make_square(
            ORANGE_MD_FROM_SOURCE("global"), {-1, -1}, {1, 1});

        UnitBuilder build;
        build.region({{outside, global}},
                     ZOrder::exterior,
                     ORANGE_MD_FROM_SOURCE("exterior"));
        build.region({{inside, hole_shape}},
                     ZOrder::hole,
                     ORANGE_MD_FROM_SOURCE("hole"));
        for (const auto* shape : {&upright, &upleft, &loleft, &loright})
        {
            build.region({{inside, *shape}},
                         ZOrder::media,
                         ORANGE_MD_FROM_SOURCE((*shape)->name()));
        }

        EXPECT_FALSE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("squares"));
    }

    VecInt check_initialization();
};

auto CoincidentSquareTest::check_initialization() -> VecInt
{
    CELER_EXPECT(this->tracker);

    std::vector<int> cells;

    // Test corners, edges, cells
    for (real_type x : {-1.0, -.5, 0.0, .5, 1.0})
    {
        for (real_type y : {-1.0, -.5, 0.0, .5, 1.0})
        {
            this->set_state({x, y, 0.0}, {1, 0, 0});
            Initialization init;
            try
            {
                init = tracker->initialize(this->state_ref());
            }
            catch (const std::exception& e)
            {
                const auto& state = this->state_ref();
                cout << "Failed to initialize at " << state.pos << " along "
                     << state.dir << ": " << e.what() << endl;
            }
            cells.push_back(init.cell ? init.cell.get() : -1);
        }
    }
    return cells;
}

//---------------------------------------------------------------------------//
//! On an upper right hole, the interior surfaces are "on" (outside)
TEST_F(CoincidentSquareTest, upper_right_coincident)
{
    this->build(
        this->make_square(ORANGE_MD_FROM_SOURCE("urhole"), {0, 0}, {1, 1}));
    CELER_ASSERT(tracker);
    EXPECT_EQ(6, tracker->num_volumes());
    EXPECT_EQ(6 + 2, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    {
        auto      result            = this->check_initialization();
        const int expected_result[] = {-1, -1, -1, -1, -1, -1, 4,  -1, 3,
                                       -1, -1, -1, -1, -1, -1, -1, 5,  -1,
                                       1,  -1, -1, -1, -1, -1, -1};
        EXPECT_VEC_EQ(expected_result, result);
    }

    {
        auto result = this->track({-2, 0.5, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "upper-left", "hole", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "lower-left.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track other direction
        auto result = this->track({2, 0.5, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "hole", "upper-left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.px", "lower-left.px", "global.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Track right along cutting plane
        auto result = this->track({-2, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "upper-left", "hole", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "lower-left.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track left along cutting plane
        auto result = this->track({2, 0, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "hole", "upper-left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.px", "lower-left.px", "global.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
        /*** END CODE ***/
    }

    {
        // Track upper right through diagonal edges
        auto result = this->track({-2, -2, 0}, {1, 1, 0});

        static const char* const expected_cells[]
            = {"exterior", "lower-left", "hole", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "lower-left.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[]
            = {1.414213562373, 1.414213562373, 1.414213562373, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

//---------------------------------------------------------------------------//
//! On a lower left hole, the interior square surfaces are outside the hole
TEST_F(CoincidentSquareTest, lower_left_coincident)
{
    this->build(
        this->make_square(ORANGE_MD_FROM_SOURCE("llhole"), {-1, -1}, {0, 0}));
    EXPECT_EQ(6, tracker->num_volumes());
    EXPECT_EQ(6 + 2, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    {
        auto      result            = this->check_initialization();
        const int expected_result[] = {-1, -1, -1, -1, -1, -1, 1,  -1, 3,
                                       -1, -1, -1, -1, -1, -1, -1, 5,  -1,
                                       2,  -1, -1, -1, -1, -1, -1};
        EXPECT_VEC_EQ(expected_result, result);
    }

    {
        // Track right through hole
        auto result = this->track({-2, -0.5, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "hole", "lower-right", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "llhole.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track left through hole
        auto result = this->track({2, -0.5, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "lower-right", "hole", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.px", "llhole.px", "global.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Track right along cutting plane
        auto result = this->track({-2, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "upper-left", "upper-right", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "llhole.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track reversed
        auto result = this->track({2, 0, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "upper-right", "upper-left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.px", "llhole.px", "global.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track upper right through diagonal edges
        auto result = this->track({-2, -2, 0}, {1, 1, 0});

        static const char* const expected_cells[]
            = {"exterior", "hole", "upper-right", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "llhole.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[]
            = {1.414213562373, 1.414213562373, 1.414213562373, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

//---------------------------------------------------------------------------//
//! On an upper right hole, the interior surfaces are "on" (outside)
TEST_F(CoincidentSquareTest, lower_right_coincident)
{
    this->build(
        this->make_square(ORANGE_MD_FROM_SOURCE("lrhole"), {0, -1}, {1, 0}));
    EXPECT_EQ(6, tracker->num_volumes());
    EXPECT_EQ(6 + 2, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    // Initialize across hole
    {
        // In lower right hole on surface
        this->set_state(
            {0.0, -0.5, 0.0}, {-1, 0, 0}, "hole", "lower-left.px", neg);
        auto init = tracker->initialize(this->state_ref());
        EXPECT_EQ(VolumeId{4}, init.cell);     // lower-left
        EXPECT_EQ(SurfaceId{6}, init.surface); // lower-left.px
    }

    {
        auto      result            = this->check_initialization();
        const int expected_result[] = {-1, -1, -1, -1, -1, -1, 4,  -1, 3,
                                       -1, -1, -1, -1, -1, -1, -1, 1,  -1,
                                       2,  -1, -1, -1, -1, -1, -1};
    }
    {
        // Track right through hole
        auto result = this->track({-2, -0.5, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "lower-left", "hole", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.mx", "lower-left.px", "global.px"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Track reversed
        auto result = this->track({2, -0.5, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"exterior", "hole", "lower-left", "exterior"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"", "global.px", "lower-left.px", "global.mx"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 1, inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Overlapping hole test.
 */
class OverlapTest : public MaskedUnitTrackerTest
{
    void SetUp()
    {
        // Overlapping spheres and enclosing sphere
        auto left = this->shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("left"), Transform{{-0.5, 0, 0}}, 1.0);
        auto right = this->shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("right"), Transform{{0.5, 0, 0}}, 1.0);
        auto parent = this->shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("parent"), Transform{{0, 0, 0}}, 3.0);

        // Two regions
        UnitBuilder build;
        build.exterior({{inside, parent}},
                       ZOrder::exterior,
                       ORANGE_MD_FROM_SOURCE("exterior"));
        build.region(
            {{inside, left}}, ZOrder::hole, ORANGE_MD_FROM_SOURCE("left"));
        build.region(
            {{inside, right}}, ZOrder::hole, ORANGE_MD_FROM_SOURCE("right"));
        build.region(
            {{inside, parent}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("inside"));
        EXPECT_FALSE(build.is_simple());
        this->make_tracker<Tracker_t>(std::move(build),
                                      ORANGE_MD_FROM_SOURCE("single surface"));
    }
};

TEST_F(OverlapTest, accessors)
{
    EXPECT_EQ(4, tracker->num_volumes());
    EXPECT_EQ(3, tracker->num_surfaces());

    // cout << celeritas::to_stream(*md, *tracker);

    std::string expected = R"rst(
======= = ======== ============================================================
Cell    Z Name     Surface logic
======= = ======== ============================================================
0       B exterior 0
                   (from ``track/test/tstMaskedUnitTracker.cc:NNNN``)
1       H left     1 ~
                   (from ``track/test/tstMaskedUnitTracker.cc:NNNN``)
2       H right    2 ~
                   (from ``track/test/tstMaskedUnitTracker.cc:NNNN``)
3       M inside   0 ~
                   (from ``track/test/tstMaskedUnitTracker.cc:NNNN``)
======= = ======== ============================================================

======= ======== ============================================================
Surface Name     Description
======= ======== ============================================================
0                Sphere: r=3
        parent.s (from ``track/test/tstMaskedUnitTracker.cc:NNNN``)
1                Sphere: r=1 at -0.5 0 0
        left.s   (from ``track/test/tstMaskedUnitTracker.cc:NNNN``)
2                Sphere: r=1 at 0.5 0 0
        right.s  (from ``track/test/tstMaskedUnitTracker.cc:NNNN``)
======= ======== ============================================================

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(OverlapTest, initialize_overlap)
{
    auto result = this->track({0, 0, 0}, {1, 1, 1});

    static const char* const expected_cells[]
        = {"Multiple overlapping cells found", "left", "right"};
    EXPECT_VEC_EQ(expected_cells, result.cells);
}

TEST_F(OverlapTest, track_overlap)
{
    {
        SCOPED_TRACE("Track rightward");
        // Enter overlapping region when exiting "left"
        auto result = this->track({-2, 0.0, 0}, {1, 0, 0});

        static const char* const expected_cells[] = {"inside",
                                                     "left",
                                                     "Multiple overlapping "
                                                     "cells found",
                                                     "right",
                                                     "left"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {"", "left.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {0.5, 1};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        SCOPED_TRACE("Track from middle");
        // Start in overlapping region
        auto result = this->track({0, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"Multiple overlapping cells found", "left", "right"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
    }
    {
        SCOPED_TRACE("Track leftward");
        // Enter overlapping region when exiting "right"
        auto result = this->track({2, 0.0, 0}, {-1, 0, 0});

        static const char* const expected_cells[] = {"inside",
                                                     "right",
                                                     "Multiple overlapping "
                                                     "cells found",
                                                     "left",
                                                     "right"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {"", "right.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' ', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {0.5, 1};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        SCOPED_TRACE("Track upward");
        auto result = this->track({0, -1.0, 0.5}, {0, 1, 0});
        static const char* const expected_cells[]
            = {"inside", "Multiple overlapping cells found", "left", "right"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {""};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const int expected_senses[] = {' '};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {0.2928932188135};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}
