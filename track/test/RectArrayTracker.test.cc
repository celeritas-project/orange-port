//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/tstRectArrayTracker.cc
 * \brief Tests for class RectArrayTracker
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../RectArrayTracker.hh"
#include "ArrayTrackerTest.hh"

#include "celeritas_test.hh"
#include "base/Face.hh"
#include "base/Future.hh"
#include "base/VectorFunctions.hh"
#include "base/Casts.hh"
#include "orange/query/RectArrayMetadata.hh"

using make_span;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class RectArrayTrackerTest
    : public orange_test::ArrayTrackerTest<celeritas::RectArrayTracker>
{
  protected:
    //// TYPE ALIASES ////
    using Grid_t         = Tracker_t::Grid_t;
    using ObjectMetadata = celeritas::ObjectMetadata;

    void make_tracker(Tracker_t::Grid_t grid, ObjectMetadata md)
    {
        CELER_EXPECT(!tracker);

        using celeritas::RectArrayMetadata;
        RectArrayMetadata::Params params;
        params.dims = grid.num_cells_dims();
        params.unit = std::move(md);
        params.bbox = {grid.low_corner(), grid.high_corner()};

        this->set_tracker(
            make_unique<Tracker_t>(std::move(grid)),
            std::make_shared<celeritas::RectArrayMetadata>(std::move(params)));
    }
};

//---------------------------------------------------------------------------//
// SINGLE CELL TESTS
class SingleCellTest : public RectArrayTrackerTest
{
    void SetUp() override
    {
        this->make_tracker({{1.0, 2.0}, {2.0, 4.0}, {3.0, 6.0}},
                           ORANGE_MD_FROM_SOURCE("single cell"));
    }
};

TEST_F(SingleCellTest, accessors)
{
    EXPECT_EQ(1, tracker->num_volumes());
    EXPECT_EQ(6, tracker->num_surfaces());

    std::string expected = R"rst(
:Dimensions: [1, 1, 1]
:x: [1, 2]
:y: [2, 4]
:z: [3, 6]

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(SingleCellTest, initialize_cell)
{
    // Initialize in center
    this->set_state({1.5, 3.0, 4.5}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize on boundaries
    this->set_state({1.0, 3.0, 4.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ("[none]", this->id_to_label(init.cell));
    EXPECT_EQ("[none]", this->id_to_label(init.surface));

    this->set_state({2.0, 3.0, 4.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ("[none]", this->id_to_label(init.cell));
    EXPECT_EQ("[none]", this->id_to_label(init.surface));

    this->set_state({2.0, 4.0, 6.0}, {1, 1, 1});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ("[none]", this->id_to_label(init.cell));
    EXPECT_EQ("[none]", this->id_to_label(init.surface));

    // Initialize outside
    this->set_state({-10, 3.0, 4.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ("[none]", this->id_to_label(init.cell));
    EXPECT_EQ("[none]", this->id_to_label(init.surface));

    this->set_state({-10, 3.0, 4.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ("[none]", this->id_to_label(init.cell));
    EXPECT_EQ("[none]", this->id_to_label(init.surface));
}

// For a single cell, we should never be able to cross into the exterior.
// (Neither surface IDs nor cells change)
TEST_F(SingleCellTest, initialize_surface)
{
    DimVector c = {0, 0, 0};
    // +X face
    this->set_state({2.0, 3.0, 4.5}, {1, 0, 0}, c, Face_t::POSX);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);

    // -Y face
    this->set_state({1.0, 2.0, 4.5}, {0, -1, 0}, c, Face_t::NEGY);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);
    EXPECT_EQ("-y", this->id_to_label(this->state_ref().surface));

    // +Z face (corner headed out)
    this->set_state({2, 4, 6}, {1, 1, 1}, c, Face_t::POSZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);
    EXPECT_EQ("+z", this->id_to_label(this->state_ref().surface));
}

// For a single cell, we should never be able to cross into the exterior.
TEST_F(SingleCellTest, intersect)
{
    const real_type eps = 1e-8;
    DimVector       c   = {0, 0, 0};
    {
        // Test from center
        static const ExpectedIntersection tests[]
            = {{{1, 0, 0}, Face_t::POSX, 0.5, c},
               {{-1, 0, 0}, Face_t::NEGX, 0.5, c},
               {{0, 1, 0}, Face_t::POSY, 1.0, c},
               {{0, -1, 0}, Face_t::NEGY, 1.0, c},
               {{0, 0, 1}, Face_t::POSZ, 1.5, c},
               {{0, 0, -1}, Face_t::NEGZ, 1.5, c}};

        this->test_intersect({1.5, 3.0, 4.5}, c, make_span(tests));
    }

    {
        // Test from off-center
        static const ExpectedIntersection tests[]
            = {{{1, 0, 0}, Face_t::POSX, 0.6, c},
               {{-1, 0, 0}, Face_t::NEGX, 0.4, c},
               {{0, 1, 0}, Face_t::POSY, 1.2, c},
               {{0, -1, 0}, Face_t::NEGY, 0.8, c},
               {{0, 0, 1}, Face_t::POSZ, 0.4, c},
               {{0, 0, -1}, Face_t::NEGZ, 2.6, c}};
        this->test_intersect({1.4, 2.8, 5.6}, c, make_span(tests));
    }

    {
        // Test off-direction from center
        static const ExpectedIntersection tests[]
            = {{{1, 1, 1}, Face_t::POSX, 0.86602540378443849, c},
               {{-1, -1, -1}, Face_t::NEGX, 0.86602540378443849, c},
               {{1, 1, 0}, Face_t::POSX, 0.70710678118654757, c},
               {{-1, -1, 0}, Face_t::NEGX, 0.70710678118654757, c},
               {{0, 1, 1}, Face_t::POSY, 1.4142135623730951, c},
               {{0, -1, -1}, Face_t::NEGY, 1.4142135623730951, c}};
        this->test_intersect({1.5, 3.0, 4.5}, c, make_span(tests));
    }
}

TEST_F(SingleCellTest, normal)
{
    DimVector c = {0, 0, 0};
    // +X face
    this->set_state({2.0, 3.0, 4.5}, {1, 0, 0}, c, Face_t::POSX);
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), tracker->normal(this->state_ref()));

    // -Y face
    this->set_state({1.0, 2.0, 4.5}, {0, -1, 0}, c, Face_t::NEGY);
    EXPECT_VEC_SOFT_EQ(Real3(0, -1, 0), tracker->normal(this->state_ref()));

    // +Z face (corner headed out)
    this->set_state({2, 4, 6}, {1, 1, 1}, c, Face_t::POSZ);
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 1), tracker->normal(this->state_ref()));
}

//---------------------------------------------------------------------------//
// MULTI-CELL TESTS
class MultiCellTest : public RectArrayTrackerTest
{
    void SetUp() override
    {
        this->make_tracker(
            {{1.0, 1.5, 2.0, 2.5}, {2.0, 4.0, 6.0, 8.0}, {3.0, 6.0, 9.0, 12.0}},
            ORANGE_MD_FROM_SOURCE("multi cell"));
    }
};

TEST_F(MultiCellTest, accessors)
{
    EXPECT_EQ(27, tracker->num_volumes());
    EXPECT_EQ(6, tracker->num_surfaces()); // Surfaces are just faces
}

TEST_F(MultiCellTest, initialize_cell)
{
    // Initialize in center
    this->set_state({1.75, 3.0, 4.5}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{1}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize way outside (reject)
    this->set_state({-15.0, 3.0, 4.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{}, init.cell);

    this->set_state({15.0, 3.0, 4.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{}, init.cell);
}

TEST_F(MultiCellTest, initialize_surface)
{
    DimVector c = {0, 0, 0};
    // Cross +X face
    this->set_state({2.0, 3.0, 4.5}, {1, 0, 0}, c, Face_t::POSX);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEGX), init.surface);
    EXPECT_EQ("{1,0,0}", this->id_to_label(init.cell));

    // Crossing -Y face leaves us in the same cell
    this->set_state({1.0, 2.0, 4.5}, {0, -1, 0}, c, Face_t::NEGY);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEGY), init.surface);

    // Crossing +Z face from corner should leave us in cell on same face
    c = {2, 2, 2};
    this->set_state({2.5, 8, 12}, {1, 1, 1}, c, Face_t::POSZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POSZ), init.surface);

    // Cross interior -X face
    c = {1, 0, 0};
    this->set_state({1.5, 3.0, 4.5}, {-1, 0, 0}, c, Face_t::NEGX);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POSX), init.surface);
}

// For a single cell, we should never be able to cross into the exterior.
TEST_F(MultiCellTest, intersect)
{
    DimVector c = {1, 1, 1};
    {
        // Test from center cell
        static const ExpectedIntersection tests[]
            = {{{1, 0, 0}, Face_t::POSX, 0.25, {2, 1, 1}},
               {{-1, 0, 0}, Face_t::NEGX, 0.25, {0, 1, 1}},
               {{0, 1, 0}, Face_t::POSY, 1.0, {1, 2, 1}},
               {{0, -1, 0}, Face_t::NEGY, 1.0, {1, 0, 1}},
               {{0, 0, 1}, Face_t::POSZ, 1.5, {1, 1, 2}},
               {{0, 0, -1}, Face_t::NEGZ, 1.5, {1, 1, 0}}};
        this->test_intersect({1.75, 5.0, 7.5}, c, make_span(tests));
    }
}

TEST_F(MultiCellTest, normal)
{
    DimVector c = {1, 0, 0};

    // +X face
    this->set_state({2.0, 3.0, 4.5}, {1, 0, 0}, c, Face_t::POSX);
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), tracker->normal(this->state_ref()));

    // -Y face
    this->set_state({1.0, 2.0, 4.5}, {0, -1, 0}, c, Face_t::NEGY);
    EXPECT_VEC_SOFT_EQ(Real3(0, -1, 0), tracker->normal(this->state_ref()));

    // +Z face (corner headed out)
    this->set_state({2, 4, 6}, {1, 1, 1}, c, Face_t::POSZ);
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 1), tracker->normal(this->state_ref()));
}
