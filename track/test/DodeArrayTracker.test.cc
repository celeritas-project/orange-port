//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/tstDodeArrayTracker.cc
 * \brief Tests for class DodeArrayTracker
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../DodeArrayTracker.hh"
#include "ArrayTrackerTest.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "base/Future.hh"
#include "base/VectorFunctions.hh"
#include "orange/query/DodeArrayMetadata.hh"

using make_span;
using constants::sqrt_two;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class DodeArrayTrackerTest
    : public orange_test::ArrayTrackerTest<celeritas::DodeArrayTracker>
{
  protected:
    //// TYPE ALIASES ////
    using ObjectMetadata = celeritas::ObjectMetadata;

    void make_tracker(real_type apothem, DimVector dims, ObjectMetadata md)
    {
        CELER_EXPECT(!tracker);

        using celeritas::DodeArrayMetadata;
        DodeArrayMetadata::Params params;
        params.apothem = apothem;
        params.dims    = dims;
        params.unit    = std::move(md);
        params.bbox    = geometria::infinite_bbox();

        this->set_tracker(
            make_unique<Tracker_t>(apothem, dims),
            std::make_shared<celeritas::DodeArrayMetadata>(std::move(params)));
    }
};

//---------------------------------------------------------------------------//
// SINGLE CELL TESTS
class SingleCellTest : public DodeArrayTrackerTest
{
    void SetUp() override
    {
        this->make_tracker(1, {1, 1, 1}, ORANGE_MD_FROM_SOURCE("single cell"));
    }
};

TEST_F(SingleCellTest, accessors)
{
    EXPECT_EQ(1, tracker->num_volumes());
    EXPECT_EQ(12, tracker->num_surfaces());

    std::string expected = R"rst(
:Type: dodecahedral array
:Apothem: 1
:Dimensions: [1, 1, 1]

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(SingleCellTest, initialize_cell)
{
    // Initialize in center
    this->set_state({0, 0, 0}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize on boundaries
    // XXX doesn't currently check for being on surface
    this->set_state({1.0, 0, 0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({-1.0, 0, 0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({0, 1.0, 0}, {1, 1, 1});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);
}

// For a single cell, we should never be able to cross into the exterior.
// (Neither surface IDs nor cells change)
// Note that `set_state` here initializes on the "outside" as though we're
// leaving the face.
TEST_F(SingleCellTest, initialize_surface)
{
    DimVector c  = {0, 0, 0};
    real_type zd = sqrt_two / 2.0;

    // +X face
    this->set_state({0, 0, 0}, {1, 0, 0}, c, Face_t::POSX);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);
    EXPECT_EQ(Sense::inside, init.sense);

    // -Y face
    this->set_state({0, 0, 0}, {0, -1, 0}, c, Face_t::NEGY);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);

    // +W face
    Real3 fw{0.5, 0.5, zd};
    Real3 fw_norm = fw;
    normalize_direction(&fw_norm);
    this->set_state(fw, fw_norm, c, Face_t::POS_A);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);

    // -W face
    this->set_state(-fw, fw_norm, c, Face_t::NEG_A);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);

    // +I face
    Real3 fi{0.5, -0.5, zd};
    Real3 fi_norm = fi;
    normalize_direction(&fi_norm);
    this->set_state(fi, fi_norm, c, Face_t::POS_D);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);
    // -I face
    this->set_state(-fi, fi_norm, c, Face_t::NEG_D);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);

    // +J face
    Real3 fj{-0.5, -0.5, zd};
    Real3 fj_norm = fj;
    normalize_direction(&fj_norm);
    this->set_state(fj, fj_norm, c, Face_t::POS_C);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);
    // -J face
    this->set_state(-fj, fj_norm, c, Face_t::NEG_C);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);

    // +K face
    Real3 fk{-0.5, 0.5, zd};
    Real3 fk_norm = fk;
    normalize_direction(&fk_norm);
    this->set_state(fk, fk_norm, c, Face_t::POS_B);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);
    // -K face
    this->set_state(-fk, fk_norm, c, Face_t::NEG_B);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);

    // Initialize on "inside" of K faces (e.g. from coincident universe)
    this->set_state(fk, fk_norm, c, Face_t::POS_B, Sense::inside);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);
    this->set_state(-fk, fk_norm, c, Face_t::NEG_B, Sense::inside);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->state_ref().surface, init.surface);
}

// For a single cell, we should never be able to cross into the exterior.
TEST_F(SingleCellTest, intersect)
{
    DimVector c = {0, 0, 0};
    {
        // Test from center
        static const ExpectedIntersection tests[] = {
            {{1, 0, 0}, Face_t::POSX, 1, c},
            {{-1, 0, 0}, Face_t::NEGX, 1, c},
            {{0, 1, 0}, Face_t::POSY, 1, c},
            {{0, -1, 0}, Face_t::NEGY, 1, c}
            // TODO add oblique face tests
        };

        this->test_intersect({0, 0, 0}, c, make_span(tests));
    }

    {
        // TODO Test from off-center
        static const ExpectedIntersection tests[] = {
            {{1, 0, 0}, Face_t::POSX, 1.4, c},
            {{-1, 0, 0}, Face_t::NEGX, 0.6, c},
            {{0, 1, 0}, Face_t::POSY, 1, c},
            {{0, -1, 0}, Face_t::NEGY, 1, c}
            // TODO add oblique face tests
        };
        this->test_intersect({-0.4, 0, 0}, c, make_span(tests));
    }

    {
        // TODO Test off-direction from center
        static const ExpectedIntersection tests[]
            = {{{1, 1, 1}, Face_t::POS_A, 1.0146118723545765, c},
               {{-1, -1, -1}, Face_t::NEG_A, 1.0146118723545765, c},
               {{-1, 1, 1}, Face_t::POS_B, 1.0146118723545765, c},
               {{1, -1, -1}, Face_t::NEG_B, 1.0146118723545765, c},
               {{1, -1, 1}, Face_t::POS_D, 1.0146118723545765, c},
               {{-1, 1, -1}, Face_t::NEG_D, 1.0146118723545765, c},
               {{-1, -1, 1}, Face_t::POS_C, 1.0146118723545765, c},
               {{1, 1, -1}, Face_t::NEG_C, 1.0146118723545765, c},
               {{1, 0.5, 0}, Face_t::POSX, 1.1180339887498949, c},
               {{-1, -0.5, 0}, Face_t::NEGX, 1.1180339887498949, c},
               {{0.5, 1, 0}, Face_t::POSY, 1.1180339887498949, c},
               {{0.5, -1, 0}, Face_t::NEGY, 1.1180339887498949, c},
               {{-1.1, -1, 0}, Face_t::NEGX, 1.3514607952107733, c}};
        this->test_intersect({0.0, 0.0, 0.0}, c, make_span(tests));
    }
}

TEST_F(SingleCellTest, normal)
{
    DimVector c  = {0, 0, 0};
    real_type zd = sqrt_two / 2.0;
    // +X face
    this->set_state({1, 0, 0}, {1, 0, 0}, c, Face_t::POSX);
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), tracker->normal(this->state_ref()));

    // -X face
    this->set_state({-1, 0, 0}, {-1, 0, 0}, c, Face_t::NEGX);
    EXPECT_VEC_SOFT_EQ(Real3(-1, 0, 0), tracker->normal(this->state_ref()));
    // -Y face
    this->set_state({0, -1.0, 0}, {0, -1, 0}, c, Face_t::NEGY);
    EXPECT_VEC_SOFT_EQ(Real3(0, -1, 0), tracker->normal(this->state_ref()));

    // +Y face
    this->set_state({0, 1, 0}, {0, 1, 0}, c, Face_t::POSY);
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), tracker->normal(this->state_ref()));

    // +W face
    Real3 fw{0.5, 0.5, zd};
    Real3 fw_norm = fw;
    normalize_direction(&fw_norm);
    this->set_state(fw, fw_norm, c, Face_t::POS_A);
    EXPECT_VEC_SOFT_EQ(fw_norm, tracker->normal(this->state_ref()));
    // -W face
    this->set_state(-fw, fw_norm, c, Face_t::NEG_A);
    EXPECT_VEC_SOFT_EQ(-fw_norm, tracker->normal(this->state_ref()));

    // +I face
    Real3 fi{0.5, -0.5, zd};
    Real3 fi_norm = fi;
    normalize_direction(&fi_norm);
    this->set_state(fi, fi_norm, c, Face_t::POS_D);
    EXPECT_VEC_SOFT_EQ(fi_norm, tracker->normal(this->state_ref()));
    // -I face
    this->set_state(-fi, fi_norm, c, Face_t::NEG_D);
    EXPECT_VEC_SOFT_EQ(-fi_norm, tracker->normal(this->state_ref()));

    // +J face
    Real3 fj{-0.5, -0.5, zd};
    Real3 fj_norm = fj;
    normalize_direction(&fj_norm);
    this->set_state(fj, fj_norm, c, Face_t::POS_C);
    EXPECT_VEC_SOFT_EQ(fj_norm, tracker->normal(this->state_ref()));
    // -J face
    this->set_state(-fj, fj_norm, c, Face_t::NEG_C);
    EXPECT_VEC_SOFT_EQ(-fj_norm, tracker->normal(this->state_ref()));

    // +K face
    Real3 fk{-0.5, 0.5, zd};
    Real3 fk_norm = fk;
    normalize_direction(&fk_norm);
    this->set_state(fk, fk_norm, c, Face_t::POS_B);
    EXPECT_VEC_SOFT_EQ(fk_norm, tracker->normal(this->state_ref()));
    // -K face
    this->set_state(-fk, fk_norm, c, Face_t::NEG_B);
    EXPECT_VEC_SOFT_EQ(-fk_norm, tracker->normal(this->state_ref()));
}

//---------------------------------------------------------------------------//
// MULTI-CELL TESTS
class MultiCellTest : public DodeArrayTrackerTest
{
  protected:
    real_type apothem_ = 2;
    real_type hheight_;
    real_type height_;

  private:
    void SetUp() override
    {
        hheight_ = sqrt_two * apothem_;
        height_  = 2 * hheight_;
        // Need some asymmetry
        this->make_tracker(
            apothem_, {5, 3, 4}, ORANGE_MD_FROM_SOURCE("multi cell"));
    }
};

TEST_F(MultiCellTest, accessors)
{
    EXPECT_EQ(60, tracker->num_volumes());
    EXPECT_EQ(12, tracker->num_surfaces()); // 12 regardless of cell count
}

TEST_F(MultiCellTest, initialize_cell)
{
    Real3 dir(1, 0, 0);
    // Initialize in center
    this->set_state({0, 0, 0}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({0.0, 2 * apothem_, 0}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 1, 0}), init.cell);

    this->set_state({2 * apothem_, 0.0, 0}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 0}), init.cell);

    // +z increments by 1
    this->set_state({apothem_, 0, apothem_}, dir);
    VolumeId zi = this->volume_id({0, 0, 1});
    init        = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);

    // Because the initial index is computed on a cuboidal [2R, 2R, 2sqrt(2)R]
    // lattice, there are 8 scenario to test.
    // 1 - the point is in the -x-y-z of initial cuboidal index (1,1,2->0,0,1
    this->set_state({1.1 * apothem_, 1.1 * apothem_, 0.1 * height_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
    // 2 - the point is in the +x-y-z of initial cuboidal index (1,0,2->0,0,1)
    this->set_state({1.1 * apothem_, 0.9 * apothem_, 0.1 * height_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
    // 3 - the point is in the -x+y-z of initial cuboidal index (0,1,2->0,0,1)
    this->set_state({0.9 * apothem_, 1.1 * apothem_, 0.1 * height_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    // 4 - the point is in the +x+y-z of initial cuboidal index (0,0,2->0,0,1)
    this->set_state({0.9 * apothem_, 0.9 * apothem_, 0.1 * height_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);

    // 5 - the point is in the -x-y+z of initial cuboidal index (1,1,1->0,0,1)
    this->set_state({1.1 * apothem_, 1.1 * apothem_, 0.9 * height_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
    // 6 - the point is in the +x-y+z of initial cuboidal index (1,0,1->0,0,1)
    this->set_state({1.1 * apothem_, 0.9 * apothem_, 0.9 * height_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
    // 7 - the point is in the -x+y+z of initial cuboidal index (0,1,1->0,0,1)
    this->set_state({0.9 * apothem_, 1.1 * apothem_, 0.9 * height_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
    // 8 - the point is in the +x+y+z of initial cuboidal index (0,0,1->0,0,1)
    this->set_state({0.9 * apothem_, 0.9 * apothem_, 0.9 * height_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);

    // test 1,1,1
    zi = this->volume_id({1, 1, 1});
    this->set_state({3 * apothem_, 3 * apothem_, hheight_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
    // test 2,2,1
    zi = this->volume_id({2, 2, 1});
    this->set_state({5 * apothem_, 5 * apothem_, hheight_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
    EXPECT_EQ("{2,2,1}", this->md->id_to_label(init.cell));
    // ensure not near boundary... really in center
    this->set_state({5.1 * apothem_, 5.1 * apothem_, hheight_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
    this->set_state({4.9 * apothem_, 4.9 * apothem_, hheight_}, dir);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(zi, init.cell);
}

TEST_F(MultiCellTest, initialize_surface)
{
    DimVector c = {0, 0, 0};
    // +X face
    this->set_state({apothem_, 0.0, 0.0}, {1, 0, 0}, c, Face_t::POSX);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEGX), init.surface);

    // +Y face
    this->set_state({0.0, apothem_, 0.0}, {0, 1, 0}, c, Face_t::POSY);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 1, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEGY), init.surface);

    // -X face
    c = {1, 0, 0};
    this->set_state({apothem_, 0.0, 0.0}, {-1, 0, 0}, c, Face_t::NEGX);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POSX), init.surface);

    // -Y face
    c = {0, 1, 0};
    this->set_state({0.0, apothem_, 0.0}, {0, -1, 0}, c, Face_t::NEGY);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POSY), init.surface);

    // +W:+X+Y+Z face
    c = {0, 0, 0};
    this->set_state({apothem_, apothem_, 0.0}, {0, 0, 1}, c, Face_t::POS_A);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_A), init.surface);

    // -W:-X-Y-Z face
    c = {0, 0, 1};
    this->set_state({apothem_, apothem_, 0.0}, {0, 0, -1}, c, Face_t::NEG_A);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_A), init.surface);
    EXPECT_EQ("+a", this->md->id_to_label(init.surface));

    // +I:+X-Y+Z face
    c = {0, 1, 0};
    this->set_state({apothem_, apothem_, 0.0}, {0, 0, 1}, c, Face_t::POS_D);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_D), init.surface);
    EXPECT_EQ("-d", this->md->id_to_label(init.surface));

    // -I:-X+Y-Z face
    c = {0, 0, 1};
    this->set_state({apothem_, apothem_, 0.0}, {0, 0, -1}, c, Face_t::NEG_D);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 1, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_D), init.surface);

    // +J:-X-Y+Z face
    c = {1, 1, 0};
    this->set_state({apothem_, apothem_, 0.0}, {0, 0, 1}, c, Face_t::POS_C);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_C), init.surface);

    // -J:+X+Y-Z face
    c = {0, 0, 1};
    this->set_state({apothem_, apothem_, 0.0}, {0, 0, -1}, c, Face_t::NEG_C);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 1, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_C), init.surface);

    // +K:-X+Y+Z face
    c = {1, 0, 0};
    this->set_state({apothem_, apothem_, 0.0}, {0, 0, 1}, c, Face_t::POS_B);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_B), init.surface);

    // -K:+X-Y-Z face
    c = {0, 0, 1};
    this->set_state({apothem_, apothem_, 0.0}, {0, 0, -1}, c, Face_t::NEG_B);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_B), init.surface);
}

TEST_F(MultiCellTest, intersect)
{
    // z direction to intersect oblique planes
    real_type zd = sqrt_two / 2.0;
    DimVector c  = {1, 1, 2};
    {
        // Test from center of even-z cell
        static const ExpectedIntersection tests[]
            //  dir         face          distance    entered cell
            = {{{1, 0, 0}, Face_t::POSX, apothem_, {2, 1, 2}},
               {{-1, 0, 0}, Face_t::NEGX, apothem_, {0, 1, 2}},
               {{0, 1, 0}, Face_t::POSY, apothem_, {1, 2, 2}},
               {{0, -1, 0}, Face_t::NEGY, apothem_, {1, 0, 2}},

               {{0.5, 0.5, zd}, Face_t::POS_A, apothem_, {1, 1, 3}},
               {{-0.5, -0.5, zd}, Face_t::POS_C, apothem_, {0, 0, 3}},
               {{-0.5, 0.5, zd}, Face_t::POS_B, apothem_, {0, 1, 3}},
               {{0.5, -0.5, zd}, Face_t::POS_D, apothem_, {1, 0, 3}},

               {{-0.5, -0.5, -zd}, Face_t::NEG_A, apothem_, {0, 0, 1}},
               {{0.5, 0.5, -zd}, Face_t::NEG_C, apothem_, {1, 1, 1}},
               {{0.5, -0.5, -zd}, Face_t::NEG_B, apothem_, {1, 0, 1}},
               {{-0.5, 0.5, -zd}, Face_t::NEG_D, apothem_, {0, 1, 1}}};
        Real3 cell_origin(4, 4, height_);
        this->test_intersect(cell_origin, c, make_span(tests));
    }
    c = {1, 1, 1};
    {
        // Test from center of odd-z cell
        static const ExpectedIntersection tests[]
            //  dir         face          distance    entered cell
            = {{{1, 0, 0}, Face_t::POSX, apothem_, {2, 1, 1}},
               {{-1, 0, 0}, Face_t::NEGX, apothem_, {0, 1, 1}},
               {{0, 1, 0}, Face_t::POSY, apothem_, {1, 2, 1}},
               {{0, -1, 0}, Face_t::NEGY, apothem_, {1, 0, 1}},

               {{0.5, 0.5, zd}, Face_t::POS_A, apothem_, {2, 2, 2}},
               {{-0.5, -0.5, zd}, Face_t::POS_C, apothem_, {1, 1, 2}},
               {{-0.5, 0.5, zd}, Face_t::POS_B, apothem_, {1, 2, 2}},
               {{0.5, -0.5, zd}, Face_t::POS_D, apothem_, {2, 1, 2}},

               {{-0.5, -0.5, -zd}, Face_t::NEG_A, apothem_, {1, 1, 0}},
               {{0.5, 0.5, -zd}, Face_t::NEG_C, apothem_, {2, 2, 0}},
               {{0.5, -0.5, -zd}, Face_t::NEG_B, apothem_, {2, 1, 0}},
               {{-0.5, 0.5, -zd}, Face_t::NEG_D, apothem_, {1, 2, 0}}};
        Real3 cell_origin(6, 6, hheight_);
        this->test_intersect(cell_origin, c, make_span(tests));
    }
}

TEST_F(MultiCellTest, normal)
{
    DimVector c = {0, 0, 1};
    // offset from coord 0,0,0 to 0,0,1
    Real3     c_offset{apothem_, apothem_, hheight_};
    real_type zd = sqrt_two / 2.0;
    // +X face
    Real3 pos{2, 0, 0};
    this->set_state(c_offset + pos, {1, 0, 0}, c, Face_t::POSX);
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), tracker->normal(this->state_ref()));

    // -X face
    pos = {-2, 0, 0};
    this->set_state(c_offset + pos, {-1, 0, 0}, c, Face_t::NEGX);
    EXPECT_VEC_SOFT_EQ(Real3(-1, 0, 0), tracker->normal(this->state_ref()));
    // -Y face
    pos = {0, -2, 0};
    this->set_state(c_offset + pos, {0, -1, 0}, c, Face_t::NEGY);
    EXPECT_VEC_SOFT_EQ(Real3(0, -1, 0), tracker->normal(this->state_ref()));

    // +Y face
    pos = {0, 2, 0};
    this->set_state(c_offset + pos, {0, 1, 0}, c, Face_t::POSY);
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), tracker->normal(this->state_ref()));

    // +W face
    Real3 fw{0.5, 0.5, zd};
    Real3 fw_norm = fw;
    normalize_direction(&fw_norm);
    this->set_state(fw + c_offset, fw_norm, c, Face_t::POS_A);
    EXPECT_VEC_SOFT_EQ(fw_norm, tracker->normal(this->state_ref()));
    // -W face
    this->set_state(-fw + c_offset, fw_norm, c, Face_t::NEG_A);
    EXPECT_VEC_SOFT_EQ(-fw_norm, tracker->normal(this->state_ref()));

    // +I face
    Real3 fi{0.5, -0.5, zd};
    Real3 fi_norm = fi;
    normalize_direction(&fi_norm);
    this->set_state(fi + c_offset, fi_norm, c, Face_t::POS_D);
    EXPECT_VEC_SOFT_EQ(fi_norm, tracker->normal(this->state_ref()));
    // -I face
    this->set_state(-fi + c_offset, fi_norm, c, Face_t::NEG_D);
    EXPECT_VEC_SOFT_EQ(-fi_norm, tracker->normal(this->state_ref()));

    // +J face
    Real3 fj{-0.5, -0.5, zd};
    Real3 fj_norm = fj;
    normalize_direction(&fj_norm);
    this->set_state(fj + c_offset, fj_norm, c, Face_t::POS_C);
    EXPECT_VEC_SOFT_EQ(fj_norm, tracker->normal(this->state_ref()));
    // -J face
    this->set_state(-fj + c_offset, fj_norm, c, Face_t::NEG_C);
    EXPECT_VEC_SOFT_EQ(-fj_norm, tracker->normal(this->state_ref()));

    // +K face
    Real3 fk{-0.5, 0.5, zd};
    Real3 fk_norm = fk;
    normalize_direction(&fk_norm);
    this->set_state(fk + c_offset, fk_norm, c, Face_t::POS_B);
    EXPECT_VEC_SOFT_EQ(fk_norm, tracker->normal(this->state_ref()));
    // -K face
    this->set_state(-fk + c_offset, fk_norm, c, Face_t::NEG_B);
    EXPECT_VEC_SOFT_EQ(-fk_norm, tracker->normal(this->state_ref()));
}
