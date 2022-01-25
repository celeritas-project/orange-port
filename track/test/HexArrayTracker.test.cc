//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/tstHexArrayTracker.cc
 * \brief Tests for class HexArrayTracker
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../HexArrayTracker.hh"
#include "ArrayTrackerTest.hh"

#include "celeritas_test.hh"
#include "base/Definitions.hh"
#include "base/Face.hh"
#include "base/Future.hh"
#include "base/VectorFunctions.hh"
#include "base/Casts.hh"
#include "orange/track/Definitions.hh"
#include "orange/query/HexArrayMetadata.hh"

using make_span;
using HexOrientation = celeritas::HexArrayTracker::Orientation;

constexpr real_type half_sqrt_three = 0.5 * constants::sqrt_three;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class HexArrayTrackerTest
    : public orange_test::ArrayTrackerTest<celeritas::HexArrayTracker>
{
    using Base = orange_test::ArrayTrackerTest<celeritas::HexArrayTracker>;

  protected:
    //// TYPE ALIASES ////
    using VecDbl         = std::vector<double>;
    using VecInt         = def::VecInt;
    using ObjectMetadata = celeritas::ObjectMetadata;

    using PlaneDimVector = Tracker_t::PlaneDimVector;
    using PlaneVector    = Tracker_t::PlaneVector;
    using Z_Edges        = VecDbl;

    void make_tracker(PlaneVector    xy_origin,
                      real_type      apothem,
                      PlaneDimVector uvdims,
                      Z_Edges        z_edges,
                      HexOrientation orientation,
                      ObjectMetadata md)
    {
        CELER_EXPECT(!tracker);

        using celeritas::HexArrayMetadata;
        HexArrayMetadata::Params params;

        params.dims = {uvdims[0], uvdims[1], z_edges.size() - 1};
        params.bbox = geometria::infinite_bbox(); // irrelevant
        params.unit = std::move(md);

        // Save arbitrary translations so we can update the unit test values
        trans_[0] = xy_origin[0];
        trans_[1] = xy_origin[1];
        trans_[2] = z_edges.front();
        for (real_type& e : z_edges)
        {
            e -= trans_[2];
        }

        this->set_tracker(
            make_unique<Tracker_t>(apothem, uvdims, z_edges, orientation),
            std::make_shared<celeritas::HexArrayMetadata>(std::move(params)));
    }

    // Test find method
    void test_find(const VecDbl& points, const VecInt& hex_coords)
    {
        CELER_EXPECT(points.size() == hex_coords.size());
        CELER_EXPECT(points.size() % 3 == 0);

        //
        auto      xyz_iter = points.begin();
        auto      uvk_iter = hex_coords.begin();
        Real3     pt;
        DimVector co;
        while (xyz_iter != points.end())
        {
            // Map data to triplets
            // Point which will be located on hex grid
            pt[0] = *xyz_iter++ - trans_[0];
            pt[1] = *xyz_iter++ - trans_[1];
            pt[2] = *xyz_iter++ - trans_[2];
            // Expected hex coordinates corresponding to this point
            co[0] = *uvk_iter++;
            co[1] = *uvk_iter++;
            co[2] = *uvk_iter++;
            //
            EXPECT_EQ(co, this->get_tracker().find(pt))
                << "Failed for location " << pt;
        }
    }

    //// Override base method to account for stupid translation ////
    // TODO: remove this and update actual test values

    // Initialize with a given pos/dir
    void set_state(const Real3& pos, const Real3& dir)
    {
        Real3 trans_pos = pos;
        trans_pos -= trans_;
        return Base::set_state(
            trans_pos, dir, VolumeId{}, SurfaceId{}, Sense::inside);
    }

    void set_state(const Real3& pos, const Real3& dir, DimVector cell)
    {
        Real3 trans_pos = pos;
        trans_pos -= trans_;
        VolumeId volume_id = this->volume_id(cell);
        return Base::set_state(
            trans_pos, dir, volume_id, SurfaceId{}, Sense::inside);
    }

    //! Override base method to account for stupid translation
    void set_state(const Real3& pos,
                   const Real3& dir,
                   DimVector    cell,
                   Face_t       face,
                   Sense        sense = Sense::outside)
    {
        CELER_EXPECT(face);
        Real3 trans_pos = pos;
        trans_pos -= trans_;
        return Base::set_state(trans_pos,
                               dir,
                               this->volume_id(cell),
                               this->surface_id(face),
                               Sense::outside);
    }

    void test_intersect(const Real3&                     pos,
                        const DimVector&                 cell,
                        span<const ExpectedIntersection> tests)
    {
        Real3 trans_pos = pos;
        trans_pos -= trans_;
        Base::test_intersect(trans_pos, cell, tests);
    }

    Real3 centroid(const DimVector& cell)
    {
        Real3 result = this->get_tracker().centroid(cell);
        result += trans_;
        return result;
    }
    // <<< End stupid overrides

  private:
    // This unit test originally supported an arbitrary "origin" for the array,
    // which was needed for GG but not ORANGE. For now, apply the transform
    // during initialization and other comparisons
    Real3 trans_;
};

//---------------------------------------------------------------------------//
// SINGLE CELL - RHOMBOIDAL FLATTOP
class SingleCellFlatTest : public HexArrayTrackerTest
{
    void SetUp() override
    {
        this->make_tracker({-0.2679491924311, -1.5},
                           1.5,
                           {1, 1},
                           {-4.0, 0.},
                           HexOrientation::flat_top,
                           ORANGE_MD_FROM_SOURCE("single cell"));
    }
};

TEST_F(SingleCellFlatTest, hex_tracker_methods)
{
    const auto* hex_tracker = &this->get_tracker();
    EXPECT_SOFT_EQ(1.5, hex_tracker->apothem());
    Z_Edges zplanes = {0, 4.0};
    EXPECT_VEC_SOFT_EQ(zplanes, hex_tracker->z_edges());
    EXPECT_EQ(HexOrientation::flat_top, hex_tracker->orientation());
    EXPECT_EQ(DimVector(1, 1, 1), hex_tracker->dims());

    // test find
    VecDbl points = {-1.09291606776104,
                     -2.67630689729359,
                     -3.60805952855753,
                     -0.90076258386996,
                     -2.65730626132298,
                     -3.64969162020404,
                     -0.45242522895066,
                     -2.98727395855945,
                     -1.33490175932897};

    VecInt expected_hex_coords = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    this->test_find(points, expected_hex_coords);

    // centroid for the given cell
    EXPECT_VEC_SOFT_EQ(Real3(-0.2679491924311, -1.5, -2.0),
                       this->centroid(DimVector{0, 0, 0}));
}

TEST_F(SingleCellFlatTest, accessors)
{
    EXPECT_EQ(1, tracker->num_volumes());
    EXPECT_EQ(8, tracker->num_surfaces());

    std::string expected = R"rst(
:Type: hexagonal array
:Orientation: flat_top
:Apothem: 1.5
:Dimensions: [1, 1, 1]
:z: [0, 4]

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(SingleCellFlatTest, initialize_cell)
{
    // Initialize in center
    this->set_state({0.0, -2.0, -2.0}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize on boundaries
    // XXX doesn't currently check for being on surface
    this->set_state({0.0, 0.0, -2.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({0.0, -3.0, -2.0}, {-1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({1.03108891196061, -0.75, -2.0}, {1, 1, 1});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize outside
    this->set_state({-4.0, -4.0, -2.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // Initialize outside
    this->set_state({10.0, 13.0, -2.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
}

// For a single cell, we should never be able to cross into the exterior.
// (Neither surface IDs nor cells change)
TEST_F(SingleCellFlatTest, initialize_surface)
{
    DimVector c = {0, 0, 0};
    // +U face
    this->set_state(
        {2.17542648054294150, -1.0, -2.0}, {1, 0, 0}, c, Face_t::POS_U);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -U face
    this->set_state({-2.0, -3.0, -2.0}, {1, 0, 0}, c, Face_t::NEG_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -U face
    this->set_state(
        {-1.7113248654051871, -2.0, -2}, {-1, 0, 0}, c, Face_t::NEG_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // +V face
    this->set_state({0.0, -3.0, -2.0}, {0, 1, 0}, c, Face_t::POS_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -V face
    this->set_state({0.0, 0.0, -2.0}, {0, -1, 0}, c, Face_t::NEG_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // +W face
    this->set_state(
        {-1.7113248654051871, -1.0, -2}, {-1, 1, 0}, c, Face_t::POS_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -W face
    this->set_state(
        {1.17542648054294150, -2.0, -2.0}, {1, -1, 0}, c, Face_t::NEG_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // +Z face
    this->set_state({0.0, 0.0, 0.0}, {0, 0, 1}, c, Face_t::POSZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -Z face
    this->set_state({0.0, 0.0, -4.0}, {0, 0, -1}, c, Face_t::NEGZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
}

TEST_F(SingleCellFlatTest, intersect)
{
    DimVector c = {0, 0, 0};
    {
        SCOPED_TRACE("from center");
        static const ExpectedIntersection tests[]
            = {{{half_sqrt_three, 0.5, 0}, Face_t::POS_U, 1.5, c},
               {{-half_sqrt_three, -0.5, 0}, Face_t::NEG_U, 1.5, c},
               {{0.0, 1.0, 0}, Face_t::POS_V, 1.5, c},
               {{0, -1, 0}, Face_t::NEG_V, 1.5, c},
               {{-half_sqrt_three, 0.5, 0}, Face_t::POS_W, 1.5, c},
               {{half_sqrt_three, -0.5, 0}, Face_t::NEG_W, 1.5, c},
               {{0.0, 0.0, 1.0}, Face_t::POSZ, 2.0, c},
               {{0.0, 0.0, -1.0}, Face_t::NEGZ, 2.0, c}};

        this->test_intersect(
            {-0.26794919243112281, -1.5, -2.0}, c, make_span(tests));
    }

    {
        SCOPED_TRACE("from off-center");
        static const ExpectedIntersection tests[] = {
            {{half_sqrt_three, 0.5, 0}, Face_t::POS_U, 1.0179491924311226, c},
            {{-half_sqrt_three, -0.5, 0}, Face_t::NEG_U, 1.982050807568877, c},
            {{0.0, 1.0, 0}, Face_t::POS_V, 1.0, c},
            {{0, -1, 0}, Face_t::NEG_V, 2.0, c},
            {{-half_sqrt_three, 0.5, 0}, Face_t::POS_W, 1.482050807568877, c},
            {{half_sqrt_three, -0.5, 0}, Face_t::NEG_W, 1.5179491924311226, c},
            {{0.0, 0.0, 1.0}, Face_t::POSZ, 1.0, c},
            {{0.0, 0.0, -1.0}, Face_t::NEGZ, 3.0, c}};

        this->test_intersect({0.0, -1.0, -1.0}, c, make_span(tests));
    }

    {
        SCOPED_TRACE("off-direction from center");
        static const ExpectedIntersection tests[]
            = {{{-1.0, -1.0, -1.0}, Face_t::NEG_U, 1.90192378864668, c},
               {{-1.0, -1.0, 0.0}, Face_t::NEG_U, 1.55291427061512, c},
               {{-1.0, -1.0, 1.0}, Face_t::NEG_U, 1.90192378864668, c},
               {{-1.0, 0.0, -1.0}, Face_t::NEG_U, 2.44948974278318, c},
               {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 1.73205080756888, c},
               {{-1.0, 0.0, 1.0}, Face_t::NEG_U, 2.44948974278318, c},
               {{-1.0, 1.0, -1.0}, Face_t::POS_W, 1.90192378864668, c},
               {{-1.0, 1.0, 0.0}, Face_t::POS_W, 1.55291427061512, c},
               {{-1.0, 1.0, 1.0}, Face_t::POS_W, 1.90192378864668, c},
               {{0.0, -1.0, -1.0}, Face_t::NEG_V, 2.12132034355964, c},
               {{0.0, -1.0, 1.0}, Face_t::NEG_V, 2.12132034355964, c},
               {{0.0, 1.0, -1.0}, Face_t::POS_V, 2.12132034355964, c},
               {{0.0, 1.0, 1.0}, Face_t::POS_V, 2.12132034355964, c},
               {{1.0, -1.0, -1.0}, Face_t::NEG_W, 1.90192378864668, c},
               {{1.0, -1.0, 0.0}, Face_t::NEG_W, 1.55291427061512, c},
               {{1.0, -1.0, 1.0}, Face_t::NEG_W, 1.90192378864668, c},
               {{1.0, 0.0, -1.0}, Face_t::POS_U, 2.44948974278318, c},
               {{1.0, 0.0, 0.0}, Face_t::POS_U, 1.73205080756888, c},
               {{1.0, 0.0, 1.0}, Face_t::POS_U, 2.44948974278318, c},
               {{1.0, 1.0, -1.0}, Face_t::POS_U, 1.90192378864668, c},
               {{1.0, 1.0, 0.0}, Face_t::POS_U, 1.55291427061512, c},
               {{1.0, 1.0, 1.0}, Face_t::POS_U, 1.90192378864668, c},
               {{0.1, 0.1, 1.0}, Face_t::POSZ, 2.01990098767242, c},
               {{0.1, 0.1, -1.0}, Face_t::NEGZ, 2.01990098767242, c}};

        this->test_intersect(
            {-0.26794919243112281, -1.5, -2.0}, c, make_span(tests));
    }
}

TEST_F(SingleCellFlatTest, normal)
{
    DimVector c = {0, 0, 0};

    // +U face
    Real3 norm = {half_sqrt_three, 0.5, 0.0};
    this->set_state({1.03108891196061, -0.75, -2.0}, norm, c, Face_t::POS_U);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // -U face
    norm = {-half_sqrt_three, -0.5, 0.0};
    this->set_state({1.5669872993927, -2.25, -2.0}, norm, c, Face_t::NEG_U);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // +V face
    norm = {0.0, 1.0, 0.0};
    this->set_state({-0.26794919371605, 0.0, -2.0}, norm, c, Face_t::POS_V);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // -V face
    norm = {0.0, -1.0, 0.0};
    this->set_state({-0.26794919371605, -3.0, -2.0}, norm, c, Face_t::NEG_V);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // +W face
    norm = {-half_sqrt_three, 0.5, 0.0};
    this->set_state({-1.56698729939271, -0.75, -2.0}, norm, c, Face_t::POS_W);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // -W face
    norm = {half_sqrt_three, -0.5, 0.0};
    this->set_state({1.03108891196061, -2.25, -2.0}, norm, c, Face_t::NEG_W);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // +Z face
    norm = {0.0, 0.0, 1.0};
    this->set_state({-1.56698729939271, -0.75, 0.0}, norm, c, Face_t::POSZ);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // -Z face
    norm = {0.0, 0.0, -1.0};
    this->set_state({1.03108891196061, -2.25, -4.0}, norm, c, Face_t::NEGZ);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));
}

//---------------------------------------------------------------------------//
// SINGLE CELL - RHOMBOIDAL POINTTOP
//---------------------------------------------------------------------------//
class SingleCellPointyTest : public HexArrayTrackerTest
{
    void SetUp() override
    {
        this->make_tracker({-0.5, -1.26794919243112},
                           1.5,
                           {1, 1},
                           {-4.0, 0.},
                           HexOrientation::pointy_top,
                           ORANGE_MD_FROM_SOURCE("single cell"));
    }
};

TEST_F(SingleCellPointyTest, hex_tracker_methods)
{
    const auto* hex_tracker = &this->get_tracker();
    EXPECT_SOFT_EQ(1.5, hex_tracker->apothem());
    EXPECT_EQ(HexOrientation::pointy_top, hex_tracker->orientation());
    EXPECT_EQ(DimVector(1, 1, 1), hex_tracker->dims());

    VecDbl points = {0.20065534520652,
                     -1.75803528388417,
                     -2.76046383858271,
                     0.65058616539844,
                     -1.95355941681844,
                     -1.77956278162247,
                     -1.67345816415685,
                     -1.48583990587848,
                     -3.89706432088658,
                     -1.48076266151469,
                     -1.66324740487441,
                     -1.80195921363156,
                     -1.03488861101803,
                     -1.66725345476091,
                     -2.48902179734416};

    VecInt expected_hex_coords = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    this->test_find(points, expected_hex_coords);

    // centroid for the given cell
    EXPECT_VEC_SOFT_EQ(Real3(-0.5, -1.2679491924311, -2.0),
                       this->centroid(DimVector(0, 0, 0)));
}

TEST_F(SingleCellPointyTest, accessors)
{
    EXPECT_EQ(1, tracker->num_volumes());
    EXPECT_EQ(8, tracker->num_surfaces());

    std::string expected = R"rst(
:Type: hexagonal array
:Orientation: pointy_top
:Apothem: 1.5
:Dimensions: [1, 1, 1]
:z: [0, 4]

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(SingleCellPointyTest, initialize_cell)
{
    // Initialize in center
    this->set_state({1, -1.2679491924311228, -2.0}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize on boundaries
    this->set_state({-2, -1.5, -2.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    this->set_state({-1.9, -1.5, -4.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    this->set_state({-1.25, -2.5669872981077808, -1.0}, {1, 1, 1});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // Initialize outside
    this->set_state({-4.0, -4.0, -2.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // Initialize outside
    this->set_state({10.0, 13.0, 4.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
}

// For a single cell, we should never be able to cross into the exterior.
// (Neither surface IDs nor cells change)
TEST_F(SingleCellPointyTest, initialize_surface)
{
    DimVector c = {0, 0, 0};
    // +U face
    this->set_state({1.0, -1.5, -2.0}, {1, 0, 0}, c, Face_t::POS_U);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -U face
    this->set_state({-2.0, -1.5, -2.0}, {0, -1, 0}, c, Face_t::NEG_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // +V face
    this->set_state(
        {0.25, 0.16506350946109644, -2.0}, {1, 1, 0}, c, Face_t::POS_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -V face
    this->set_state(
        {-1.25, -2.5669872981077808, -1.0}, {-1, -1, 0}, c, Face_t::NEG_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // +W face
    this->set_state(
        {-1.25, 0.16506350946109644, -2.0}, {-1, 1, 0}, c, Face_t::POS_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -W face
    this->set_state(
        {0.25, -2.5669872981077808, -1.0}, {1, -1, 0}, c, Face_t::NEG_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // +Z face
    this->set_state({-1.0, -1.5, 0.0}, {1, 0, 1}, c, Face_t::POSZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);

    // -Z face
    this->set_state({-1.0, -1.5, -4.0}, {1, 0, -1}, c, Face_t::NEGZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(VolumeId{0}, init.cell);
}

TEST_F(SingleCellPointyTest, intersect)
{
    DimVector c = {0, 0, 0};
    {
        SCOPED_TRACE("from center");
        static const ExpectedIntersection tests[]
            = {{{1.0, 0.0, 0.0}, Face_t::POS_U, 1.5, c},
               {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 1.5, c},
               {{0.5, half_sqrt_three, 0.0}, Face_t::POS_V, 1.5, c},
               {{-0.5, -half_sqrt_three, 0.0}, Face_t::NEG_V, 1.5, c},
               {{-0.5, half_sqrt_three, 0.0}, Face_t::POS_W, 1.5, c},
               {{0.5, -half_sqrt_three, 0.0}, Face_t::NEG_W, 1.5, c},
               {{0.0, 0.0, 1.0}, Face_t::POSZ, 2.0, c},
               {{0.0, 0.0, -1.0}, Face_t::NEGZ, 2.0, c}};

        this->test_intersect(
            {-0.5, -1.2679491924311, -2.0}, c, make_span(tests));
    }

    {
        SCOPED_TRACE("from off-center");
        static const ExpectedIntersection tests[] = {
            {{1.0, 0.0, 0.0}, Face_t::POS_U, 1.0, c},
            {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 2.0, c},
            {{0.5, half_sqrt_three, 0.0}, Face_t::POS_V, 1.0179491924311228, c},
            {{-0.5, -half_sqrt_three, 0.0}, Face_t::NEG_V, 1.9820508075688772, c},
            {{-0.5, half_sqrt_three, 0.0}, Face_t::POS_W, 1.5179491924311221, c},
            {{0.5, -half_sqrt_three, 0.0}, Face_t::NEG_W, 1.4820508075688779, c},
            {{0.0, 0.0, 1.0}, Face_t::POSZ, 1.0, c},
            {{0.0, 0.0, -1.0}, Face_t::NEGZ, 3.0, c}};

        this->test_intersect({0.0, -1.0, -1.0}, c, make_span(tests));
    }

    {
        SCOPED_TRACE("off-direction from center");
        static const ExpectedIntersection tests[]
            = {{{-1.0, -1.0, -1.0}, Face_t::NEG_V, 1.90192378864668, c},
               {{-1.0, -1.0, 0.0}, Face_t::NEG_V, 1.55291427061512, c},
               {{-1.0, -1.0, 1.0}, Face_t::NEG_V, 1.90192378864668, c},
               {{-1.0, 0.0, -1.0}, Face_t::NEG_U, 2.12132034355964, c},
               {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 1.50000000000000, c},
               {{-1.0, 0.0, 1.0}, Face_t::NEG_U, 2.12132034355964, c},
               {{-1.0, 1.0, -1.0}, Face_t::POS_W, 1.90192378864668, c},
               {{-1.0, 1.0, 0.0}, Face_t::POS_W, 1.55291427061512, c},
               {{-1.0, 1.0, 1.0}, Face_t::POS_W, 1.90192378864668, c},
               {{0.0, -1.0, -1.0}, Face_t::NEG_W, 2.44948974278318, c},
               {{0.0, -1.0, 1.0}, Face_t::NEG_W, 2.44948974278318, c},
               {{0.0, 1.0, -1.0}, Face_t::POS_W, 2.44948974278318, c},
               {{0.0, 1.0, 1.0}, Face_t::POS_W, 2.44948974278318, c},
               {{1.0, -1.0, -1.0}, Face_t::NEG_W, 1.90192378864668, c},
               {{1.0, -1.0, 0.0}, Face_t::NEG_W, 1.55291427061512, c},
               {{1.0, -1.0, 1.0}, Face_t::NEG_W, 1.90192378864668, c},
               {{1.0, 0.0, -1.0}, Face_t::POS_U, 2.12132034355964, c},
               {{1.0, 0.0, 0.0}, Face_t::POS_U, 1.50000000000000, c},
               {{1.0, 0.0, 1.0}, Face_t::POS_U, 2.12132034355964, c},
               {{1.0, 1.0, -1.0}, Face_t::POS_V, 1.90192378864668, c},
               {{1.0, 1.0, 0.0}, Face_t::POS_V, 1.55291427061512, c},
               {{1.0, 1.0, 1.0}, Face_t::POS_V, 1.90192378864668, c},
               {{0.1, 0.1, 1.0}, Face_t::POSZ, 2.01990098767242, c},
               {{0.1, 0.1, -1.0}, Face_t::NEGZ, 2.01990098767242, c}};

        this->test_intersect(
            {-0.5, -1.2679491924311, -2.0}, c, make_span(tests));
    }
}

TEST_F(SingleCellPointyTest, normal)
{
    DimVector c = {0, 0, 0};

    // +U face
    Real3 norm = {1.0, 0.0, 0.0};
    this->set_state({1.0, -1.26794922351837, -2.0}, norm, c, Face_t::POS_U);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // -U face
    norm = {-1.0, 0.0, 0.0};
    this->set_state({-2.0, -1.26794922351837, -2.0}, norm, c, Face_t::NEG_U);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // +V face
    norm = {0.5, half_sqrt_three, 0.0};
    this->set_state({0.25, 0.03108888215829, -2.0}, norm, c, Face_t::POS_V);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // -V face
    norm = {-0.5, -half_sqrt_three, 0.0};
    this->set_state({-1.25, -2.56698732919503, -2.0}, norm, c, Face_t::NEG_V);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // +W face
    norm = {-0.5, half_sqrt_three, 0.0};
    this->set_state({-1.25, 0.03108888215829, -2.0}, norm, c, Face_t::POS_W);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // -W face
    norm = {0.5, -half_sqrt_three, 0.0};
    this->set_state({0.25, -2.56698732919503, -2.0}, norm, c, Face_t::NEG_W);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // +Z face
    norm = {0.0, 0.0, 1.0};
    this->set_state({0.0, 0.0, 0.0}, norm, c, Face_t::POSZ);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));

    // -Z face
    norm = {0.0, 0.0, -1.0};
    this->set_state({0.0, 0.0, -4.0}, norm, c, Face_t::NEGZ);
    EXPECT_VEC_SOFT_EQ(norm, tracker->normal(this->state_ref()));
}

//---------------------------------------------------------------------------//
// MULTI CELL - RHOMBOIDAL FLATTOP
class MultiCellFlatTest : public HexArrayTrackerTest
{
    void SetUp() override
    {
        this->make_tracker({0.309401076758503, -1.0},
                           2.0,
                           {4, 2},
                           {-1., 0., 4.0}, // Z edges
                           HexOrientation::flat_top,
                           ORANGE_MD_FROM_SOURCE("multi cell"));
    }
};

TEST_F(MultiCellFlatTest, hex_tracker_methods)
{
    const auto* hex_tracker = &this->get_tracker();
    EXPECT_SOFT_EQ(2.0, hex_tracker->apothem());
    EXPECT_SOFT_EQ(2.3094010767585029, hex_tracker->circumradius());
    EXPECT_VEC_SOFT_EQ(Z_Edges({0.0, 1.0, 5.0}), hex_tracker->z_edges());

    EXPECT_EQ(DimVector(4, 2, 2), hex_tracker->dims());
    EXPECT_EQ(2 * 4 + 1 * 2 + 0, this->volume_id({2, 1, 0}).unchecked_get());
    EXPECT_EQ("{0,1,0}", this->md->id_to_label(this->volume_id({0, 1, 0})));

    VecDbl points = {-0.83258839819136, -2.11275132208602, -0.21389452438777,
                     0.66443386212080,  -1.77905265089150, -0.96152521359746,
                     -0.81119988052630, 0.65012278115803,  -0.23611155617081,
                     1.03926312188856,  -0.26823054624334, -0.41304375711174,
                     2.60270221427120,  1.01297011253604,  -0.58204205317372,
                     -1.05339317456804, 3.38808145171193,  -0.12922177108493,
                     0.21864012613663,  2.81976624920128,  -0.75822424022966,
                     4.15708079622111,  2.39618069366394,  -0.71208391509570,
                     5.95247545966125,  2.84548004338293,  -0.12577826066558,
                     8.14340313818122,  3.31465172507787,  -0.17705370612380,
                     2.56276772970112,  5.12979631247900,  -0.73014538262607,
                     5.98726019323569,  4.22787191603037,  -0.54211070816543,
                     7.02218759597628,  4.05175498473891,  -0.32804349522298,
                     10.03843917910891, 4.80557772415691,  -0.51260358417207,
                     12.32018260778825, 4.57116732429505,  -0.28301289855972,
                     2.96100319107898,  6.98315162451742,  -0.59339361575016,
                     4.56112601377016,  6.80836091795243,  -0.38475768609805,
                     6.79592032385689,  7.95618608119722,  -0.57211230677586,
                     9.30809826115918,  6.81960767186729,  -0.00014213911693,
                     11.78922961573489, 6.35299533348175,  -0.95903833579065,
                     9.76247916427770,  9.74020787451407,  -0.37480513265161,
                     11.59365664866402, 9.52459053364930,  -0.09796725392416};

    VecInt expected_hex_coords
        = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1,
           0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0,
           0, 1, 1, 0, 1, 1, 0, 2, 1, 0, 2, 1, 0, 3, 0, 0, 3, 1, 0, 3, 1, 0};

    this->test_find(points, expected_hex_coords);

    // centroid for the given cell
    auto r = hex_tracker->circumradius();
    EXPECT_VEC_SOFT_EQ(Real3(4 * r - 2.0, 7.0, -0.5),
                       this->centroid(DimVector{2, 1, 0}));
}

TEST_F(MultiCellFlatTest, initialize_cell)
{
    // Initialize in center
    this->set_state({0.0, 0.0, -0.5}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({3.5, 1.0, -0.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 0}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({7.5, 7.0, -0.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({2, 1, 0}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);
    EXPECT_EQ("{2,1,0}", this->md->id_to_label(init.cell));

    this->set_state({0.0, 0.0, 2.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize on outer boundaries
    this->set_state({0.0, -3.0, -1.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);

    this->set_state({0.0, -3.0, 4.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);

    this->set_state({10.0, 3.0, 4.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({3, 0, 1}), init.cell);

    // Initialize outside
    this->set_state({-4.0, -3.0, -1.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);

    this->set_state({-4.0, 3.0, 4.0}, {1, 1, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 1, 1}), init.cell);

    this->set_state({16.0, 3.5, 4.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({3, 0, 1}), init.cell);

    this->set_state({-1.6, 1.7, -0.1}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 1, 0}), init.cell);
}

// For a single cell, we should never be able to cross into the exterior.
// (Neither surface IDs nor cells change)
TEST_F(MultiCellFlatTest, initialize_surface)
{
    // +U face (inside the array)
    DimVector c = {0, 0, 0};
    this->set_state(
        {2.0414518843273806, 0.0, -0.5}, {1, 0, 0}, c, Face_t::POS_U);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_U), init.surface);
    EXPECT_EQ("-u", this->md->id_to_label(init.surface));

    c = {3, 0, 0};
    // +U face (headed outside)
    this->set_state(
        {12.433756729740645, 6.0, -0.5}, {1, 0, 0}, c, Face_t::POS_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // -U face (inside the array)
    c = {1, 0, 1};
    this->set_state(
        {2.0414518843273806, 0.0, 0.5}, {-1, 0, 0}, c, Face_t::NEG_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_U), init.surface);

    c = {0, 0, 1};
    // -U face (headed outside)
    this->set_state(
        {-1.4226497308103743, -2.0, 0.5}, {1, 0, 0}, c, Face_t::NEG_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // +V face (inside the array)
    c = {3, 0, 1};
    this->set_state({11.0, 7.0, 0.5}, {1, 1, 0}, c, Face_t::POS_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({3, 1, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_V), init.surface);

    c = {3, 1, 1};
    // +V face (headed outside)
    this->set_state({11.0, 11.0, 0.5}, {1, 1, 0}, c, Face_t::POS_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // -V face (inside the array)
    c = {2, 1, 0};
    this->set_state({7.0, 5.0, 0.5}, {0, -1, 0}, c, Face_t::NEG_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({2, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_V), init.surface);

    c = {1, 0, 1};
    // -V face (headed outside)
    this->set_state({4.0, -1.0, 0.5}, {1, 0, 0}, c, Face_t::NEG_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // +W face (inside the array)
    c = {3, 0, 1};
    this->set_state(
        {8.9696551146028902, 6.0, 0.5}, {-1, 1, 0}, c, Face_t::POS_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({2, 1, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_W), init.surface);

    c = {1, 1, 0};
    // +W face (headed outside)
    this->set_state(
        {2.0414518843273806, 6.0, -0.5}, {-1, 1, 0}, c, Face_t::POS_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // -W face (inside the array)
    c = {0, 1, 1};
    this->set_state(
        {2.0414518843273806, 2.0, 0.5}, {1, -1, 0}, c, Face_t::NEG_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_W), init.surface);
    EXPECT_EQ("+w", this->md->id_to_label(init.surface));

    c = {2, 0, 0};
    // -W face (headed outside)
    this->set_state(
        {8.9696551146028902, 2.0, -0.5}, {1, -1, 0}, c, Face_t::NEG_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // +Z face
    c = {0, 0, 0};
    this->set_state({0.0, 0.0, 0.0}, {0, 0, 1}, c, Face_t::POSZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEGZ), init.surface);

    // +Z face (headed outside)
    c = {0, 0, 1};
    this->set_state({0.0, 0.0, 4.8}, {0, 0, 1}, c, Face_t::POSZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // -Z face
    c = {0, 0, 1};
    this->set_state({0.0, 0.0, 0.0}, {0, 0, -1}, c, Face_t::NEGZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POSZ), init.surface);

    // -Z face (headed outside)
    c = {0, 0, 0};
    this->set_state({0.0, 0.0, -1.0}, {0, 0, -1}, c, Face_t::NEGZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);
}

TEST_F(MultiCellFlatTest, intersect)
{
    const auto* hex_tracker = &this->get_tracker();
    DimVector   c           = {1, 1, 1};
    {
        SCOPED_TRACE("from center");
        static const ExpectedIntersection tests[]
            = {{{half_sqrt_three, 0.5, 0}, Face_t::POS_U, 2.0, {2, 1, 1}},
               {{-half_sqrt_three, -0.5, 0}, Face_t::NEG_U, 2.0, {0, 1, 1}},
               {{0.0, 1.0, 0}, Face_t::POS_V, 2.0, c},
               {{0, -1, 0}, Face_t::NEG_V, 2.0, {1, 0, 1}},
               {{-half_sqrt_three, 0.5, 0}, Face_t::POS_W, 2.0, c},
               {{half_sqrt_three, -0.5, 0}, Face_t::NEG_W, 2.0, {2, 0, 1}},
               {{0.0, 0.0, 1.0}, Face_t::POSZ, 2.0, c},
               {{0.0, 0.0, -1.0}, Face_t::NEGZ, 2.0, {1, 1, 0}}};

        this->test_intersect(this->centroid(c), c, make_span(tests));
    }

    c = {2, 1, 0};
    {
        SCOPED_TRACE("from off-center");
        static const ExpectedIntersection tests[] = {
            {{half_sqrt_three, 0.5, 0},
             Face_t::POS_U,
             2.4820545354643317,
             {3, 1, 0}},
            {{-half_sqrt_three, -0.5, 0},
             Face_t::NEG_U,
             1.5179454645356674,
             {1, 1, 0}},
            {{0.0, 1.0, 0}, Face_t::POS_V, 2.5, c},
            {{0, -1, 0}, Face_t::NEG_V, 1.5, {2, 0, 0}},
            {{-half_sqrt_three, 0.5, 0}, Face_t::POS_W, 2.0179454645356674, c},
            {{half_sqrt_three, -0.5, 0},
             Face_t::NEG_W,
             1.9820545354643315,
             {3, 0, 0}},
            {{0.0, 0.0, 1.0}, Face_t::POSZ, 0.7, {2, 1, 1}},
            {{0.0, 0.0, -1.0}, Face_t::NEGZ, 0.3, c}};

        this->test_intersect({6.96965081, 6.5, -0.7}, c, make_span(tests));
    }

    c = {2, 1, 0};
    {
        SCOPED_TRACE("off-direction from center");
        static const ExpectedIntersection tests[]
            = {{{-1.0, -1.0, -1.0}, Face_t::NEGZ, 0.86602540378444, c},
               {{-1.0, -1.0, 0.0}, Face_t::NEG_U, 2.07055236082017, {1, 1, 0}},
               {{-1.0, -1.0, 1.0}, Face_t::POSZ, 0.86602540378444, {2, 1, 1}},
               {{-1.0, 0.0, -1.0}, Face_t::NEGZ, 0.70710678118655, c},
               {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 2.30940107675850, {1, 1, 0}},
               {{-1.0, 0.0, 1.0}, Face_t::POSZ, 0.70710678118655, {2, 1, 1}},
               {{-1.0, 1.0, -1.0}, Face_t::NEGZ, 0.86602540378444, c},
               {{-1.0, 1.0, 0.0}, Face_t::POS_W, 2.07055236082017, c},
               {{-1.0, 1.0, 1.0}, Face_t::POSZ, 0.86602540378444, {2, 1, 1}},
               {{0.0, -1.0, -1.0}, Face_t::NEGZ, 0.70710678118655, c},
               {{0.0, -1.0, 1.0}, Face_t::POSZ, 0.70710678118655, {2, 1, 1}},
               {{0.0, 1.0, -1.0}, Face_t::NEGZ, 0.70710678118655, c},
               {{0.0, 1.0, 1.0}, Face_t::POSZ, 0.70710678118655, {2, 1, 1}},
               {{1.0, -1.0, -1.0}, Face_t::NEGZ, 0.86602540378444, c},
               {{1.0, -1.0, 0.0}, Face_t::NEG_W, 2.07055236082017, {3, 0, 0}},
               {{1.0, -1.0, 1.0}, Face_t::POSZ, 0.86602540378444, {2, 1, 1}},
               {{1.0, 0.0, -1.0}, Face_t::NEGZ, 0.70710678118655, c},
               {{1.0, 0.0, 0.0}, Face_t::POS_U, 2.30940107675850, {3, 1, 0}},
               {{1.0, 0.0, 1.0}, Face_t::POSZ, 0.70710678118655, {2, 1, 1}},
               {{1.0, 1.0, -1.0}, Face_t::NEGZ, 0.86602540378444, c},
               {{1.0, 1.0, 0.0}, Face_t::POS_U, 2.07055236082017, {3, 1, 0}},
               {{1.0, 1.0, 1.0}, Face_t::POSZ, 0.86602540378444, {2, 1, 1}},
               {{0.1, 0.1, 1.0}, Face_t::POSZ, 0.50497524691810, {2, 1, 1}},
               {{0.1, 0.1, -1.0}, Face_t::NEGZ, 0.50497524691810, c}};

        this->test_intersect(this->centroid(c), c, make_span(tests));
    }
}

//---------------------------------------------------------------------------//
// MULTI CELL - RHOMBOIDAL POINTYTOP
class MultiCellPointyTest : public HexArrayTrackerTest
{
    void SetUp() override
    {
        this->make_tracker({10.0, 11.5470053837925},
                           10.0,
                           {3, 5},
                           {-2., 0., 2.0, 4.0, 6.0}, // Z edges
                           HexOrientation::pointy_top,
                           ORANGE_MD_FROM_SOURCE("multi cell"));
    }
};

TEST_F(MultiCellPointyTest, hex_tracker_methods)
{
    const auto* hex_tracker = &this->get_tracker();
    EXPECT_SOFT_EQ(10.0, hex_tracker->apothem());
    EXPECT_SOFT_EQ(11.547005383792515, hex_tracker->circumradius());

    EXPECT_EQ(DimVector(3, 5, 4), hex_tracker->dims());

    VecDbl points = {2.63448124935602,  7.40927267358990,  1.20927525219324,
                     20.38802922386129, 9.09345410321055,  3.73183108169887,
                     37.28143337543243, 5.45373712825263,  3.84610554545089,
                     52.59488775402477, 3.12686810216914,  2.52786383294713,
                     59.55444570590302, 11.07106959228737, -0.20850725691541,
                     4.24562245096277,  18.83966617564059, -1.23727787546493,
                     21.17920180027001, 18.25534597294517, -1.00796482519734,
                     39.57775116966663, 18.56276453486939, -1.20300058228481,
                     51.96233648784697, 15.47120809027053, 2.73849531142861,
                     59.82434162848259, 25.69658575401204, 1.14445058072085,
                     25.61211152747560, 43.30027678589562, -1.39122216023090,
                     39.03133919240241, 34.65650153269046, 2.47828093155178,
                     51.73358297109323, 39.60382267555526, 1.77758634668854,
                     64.15971032057861, 41.43743752826455, 0.20370872763062,
                     24.17886868181488, 52.63883128896106, -1.32339796793158,
                     30.43692606270208, 58.43722242416953, -1.71008444000344,
                     52.11758058735440, 50.05680377334804, -1.66321562262081,
                     66.92608886812936, 56.05586118490811, 2.49451892434142,
                     73.47498664627408, 46.88511100606627, -0.41356920590029,
                     89.91487952387601, 58.82021262081787, 1.01590539772798,
                     37.74880262685700, 67.07043974372527, 1.93445441010837,
                     44.25406464826602, 70.43851377642093, 2.54977883777837,
                     59.20845089580928, 62.94370211722446, 3.13534855069321,
                     75.30212193300046, 63.22498461018266, -1.86473690466777,
                     92.65315202171551, 71.19143240133197, -0.11242795152145,
                     46.14637670961022, 82.05556137840578, 3.46664347016029,
                     59.40205716261166, 80.41191746136340, -0.13496840874308,
                     76.60675013617782, 78.11997668809437, 1.48358460384381,
                     94.30133305744465, 87.85724093409748, 0.63842244909591};

    VecInt expected_hex_coords
        = {0, 0, 1, 1, 0, 2, 1, 0, 2, 2, 0, 2, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1,
           1, 0, 2, 0, 2, 2, 1, 1, 0, 2, 0, 1, 1, 2, 1, 2, 1, 2, 2, 1, 0, 2,
           0, 0, 3, 0, 1, 2, 0, 1, 3, 2, 2, 2, 0, 2, 3, 1, 0, 3, 1, 0, 3, 2,
           1, 3, 2, 2, 3, 0, 2, 4, 0, 0, 4, 2, 0, 4, 0, 1, 4, 1, 2, 4, 1};

    this->test_find(points, expected_hex_coords);
}

TEST_F(MultiCellPointyTest, accessors)
{
    EXPECT_EQ(60, tracker->num_volumes());
    EXPECT_EQ(8, tracker->num_surfaces());

    std::string expected = R"rst(
:Type: hexagonal array
:Orientation: pointy_top
:Apothem: 10
:Dimensions: [3, 5, 4]
:z: [0, 2, 4, 6, 8]

)rst";
    EXPECT_EQ(expected, this->describe_md());
}

TEST_F(MultiCellPointyTest, initialize_cell)
{
    // Initialize in center
    this->set_state({10.0, 10.0, 0.0}, {1, 0, 0});
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({50.0, 10.0, -0.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({2, 0, 0}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({60.0, 60.0, 5.5}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 3, 3}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    this->set_state({50.0, 80.0, 5.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 4, 3}), init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);

    // Initialize on boundaries
    this->set_state({0.0, 10.0, -1.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);

    this->set_state({60.0, 10.0, 4.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({2, 0, 3}), init.cell);

    this->set_state({40.0, 80.0, 6.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 4, 3}), init.cell);

    // Initialize outside
    this->set_state({-1.0, -1.0, -3.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 0}), init.cell);

    this->set_state({65.0, 3.0, 4.0}, {1, 1, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({2, 0, 3}), init.cell);

    this->set_state({80.0, 100.5, 3.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 4, 2}), init.cell);

    this->set_state({10.0, 10.0, 25.0}, {1, 0, 0});
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 3}), init.cell);
}

// For a single cell, we should never be able to cross into the exterior.
// (Neither surface IDs nor cells change)
TEST_F(MultiCellPointyTest, initialize_surface)
{
    // +U face (inside the array)
    DimVector c = {0, 0, 0};
    this->set_state({20.0, 10.0, -0.5}, {1, 0, 0}, c, Face_t::POS_U);
    auto init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 0}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_U), init.surface);

    c = {2, 0, 0};
    // +U face (headed outside)
    this->set_state({80.0, 50.0, -0.5}, {1, 0, 0}, c, Face_t::POS_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // -U face (inside the array)
    c = {1, 2, 1};
    this->set_state({40.0, 50.0, 0.5}, {-1, 0, 0}, c, Face_t::NEG_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 2, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_U), init.surface);

    c = {0, 4, 1};
    // -U face (headed outside)
    this->set_state({40.0, 80.0, 0.5}, {1, 0, 0}, c, Face_t::NEG_U);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // +V face (inside the array)
    c = {1, 2, 2};
    this->set_state(
        {55.0, 54.848275573014455, 0.5}, {1, 1, 0}, c, Face_t::POS_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 3, 2}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_V), init.surface);

    c = {1, 4, 2};
    // +V face (headed outside)
    this->set_state(
        {75.0, 89.48929172439199, 0.5}, {1, 1, 0}, c, Face_t::POS_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // -V face (inside the array)
    c = {0, 3, 3};
    this->set_state(
        {35.0, 54.848275573014455, 3.0}, {0, -1, 0}, c, Face_t::NEG_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 2, 3}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_V), init.surface);

    c = {2, 0, 3};
    // -V face (headed outside)
    this->set_state(
        {45.0, 2.8867513459481282, 4.5}, {1, 0, 0}, c, Face_t::NEG_V);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // +W face (inside the array)
    c = {2, 0, 1};
    this->set_state(
        {45.0, 20.207259421636898, 0.5}, {-1, 1, 0}, c, Face_t::POS_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 1, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEG_W), init.surface);

    c = {0, 0, 0};
    // +W face (headed outside)
    this->set_state(
        {5.0, 20.207259421636898, -0.5}, {-1, 1, 0}, c, Face_t::POS_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // -W face (inside the array)
    c = {0, 1, 3};
    this->set_state(
        {25.0, 20.207259421636898, 4.5}, {1, -1, 0}, c, Face_t::NEG_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({1, 0, 3}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POS_W), init.surface);

    c = {2, 1, 2};
    // -W face (headed outside)
    this->set_state(
        {65.0, 20.207259421636898, 3.1}, {1, -1, 0}, c, Face_t::NEG_W);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // +Z face
    c = {0, 0, 1};
    this->set_state({10.0, 10.0, 2.0}, {0, 0, 1}, c, Face_t::POSZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 2}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::NEGZ), init.surface);

    // +Z face (headed outside)
    c = {0, 0, 3};
    this->set_state({10.0, 10.0, 6.0}, {0, 0, 1}, c, Face_t::POSZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);

    // -Z face
    c = {0, 0, 2};
    this->set_state({0.0, 0.0, 2.0}, {0, 0, -1}, c, Face_t::NEGZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id({0, 0, 1}), init.cell);
    EXPECT_EQ(this->surface_id(Face_t::POSZ), init.surface);

    // -Z face (headed outside)
    c = {0, 0, 0};
    this->set_state({10.0, 10.0, -2.0}, {0, 0, -1}, c, Face_t::NEGZ);
    init = tracker->initialize(this->state_ref());
    EXPECT_EQ(this->volume_id(c), init.cell);
}

TEST_F(MultiCellPointyTest, intersect)
{
    const auto* hex_tracker = &this->get_tracker();

    DimVector c = {0, 0, 0};
    {
        SCOPED_TRACE("from center, exterior cell");
        static const ExpectedIntersection tests[]
            = {{{1.0, 0.0, 0.0}, Face_t::POS_U, 10.0, {1, 0, 0}},
               {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 10.0, c},
               {{0.5, half_sqrt_three, 0.0}, Face_t::POS_V, 10.0, {0, 1, 0}},
               {{-0.5, -half_sqrt_three, 0.0}, Face_t::NEG_V, 10.0, c},
               {{-0.5, half_sqrt_three, 0.0}, Face_t::POS_W, 10.0, c},
               {{0.5, -half_sqrt_three, 0.0}, Face_t::NEG_W, 10.0, c},
               {{0.0, 0.0, 1.0}, Face_t::POSZ, 1.0, {0, 0, 1}},
               {{0.0, 0.0, -1.0}, Face_t::NEGZ, 1.0, c}};

        this->test_intersect(
            {10.0, 11.547005383792515, -1.0}, c, make_span(tests));
    }

    c = {1, 2, 2};
    {
        SCOPED_TRACE("from center, interior cell");
        static const ExpectedIntersection tests[]
            = {{{1.0, 0.0, 0.0}, Face_t::POS_U, 10.0, {2, 2, 2}},
               {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 10.0, {0, 2, 2}},
               {{0.5, half_sqrt_three, 0.0}, Face_t::POS_V, 10.0, {1, 3, 2}},
               {{-0.5, -half_sqrt_three, 0.0}, Face_t::NEG_V, 10.0, {1, 1, 2}},
               {{-0.5, half_sqrt_three, 0.0}, Face_t::POS_W, 10.0, {0, 3, 2}},
               {{0.5, -half_sqrt_three, 0.0}, Face_t::NEG_W, 10.0, {2, 1, 2}},
               {{0.0, 0.0, 1.0}, Face_t::POSZ, 1.0, {1, 2, 3}},
               {{0.0, 0.0, -1.0}, Face_t::NEGZ, 1.0, {1, 2, 1}}};

        this->test_intersect(this->centroid(c), c, make_span(tests));
    }

    c = {0, 0, 0};
    {
        SCOPED_TRACE("from off-center, exterior cell");
        static const ExpectedIntersection tests[] = {
            {{1.0, 0.0, 0.0}, Face_t::POS_U, 9.00000000000000, {1, 0, 0}},
            {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 11.00000000000000, c},
            {{0.5, half_sqrt_three, 0.0},
             Face_t::POS_V,
             8.63397459621556,
             {0, 1, 0}},
            {{-0.5, -half_sqrt_three, 0.0}, Face_t::NEG_V, 11.36602540378444, c},
            {{-0.5, half_sqrt_three, 0.0}, Face_t::POS_W, 9.63397459621556, c},
            {{0.5, -half_sqrt_three, 0.0}, Face_t::NEG_W, 10.36602540378444, c},
            {{0.0, 0.0, 1.0}, Face_t::POSZ, 0.35, {0, 0, 1}},
            {{0.0, 0.0, -1.0}, Face_t::NEGZ, 1.65, c}};

        this->test_intersect(
            {11.0, 12.547005383792515, -0.35}, c, make_span(tests));
    }

    c = {2, 3, 3};
    {
        SCOPED_TRACE("from off-direction from center, exterior cell");
        static const ExpectedIntersection tests[] = {
            {{-1.0, -1.0, -1.0}, Face_t::NEGZ, 1.73205080756888, {2, 3, 2}},
            {{-1.0, -1.0, 0.0}, Face_t::NEG_V, 10.35276180410083, {2, 2, 3}},
            {{-1.0, -1.0, 1.0}, Face_t::POSZ, 1.73205080756888, c},
            {{-1.0, 0.0, -1.0}, Face_t::NEGZ, 1.41421356237310, {2, 3, 2}},
            {{-1.0, 0.0, 0.0}, Face_t::NEG_U, 10.00000000000000, {1, 3, 3}},
            {{-1.0, 0.0, 1.0}, Face_t::POSZ, 1.41421356237310, c},
            {{-1.0, 1.0, -1.0}, Face_t::NEGZ, 1.73205080756888, {2, 3, 2}},
            {{-1.0, 1.0, 0.0}, Face_t::POS_W, 10.35276180410083, {1, 4, 3}},
            {{-1.0, 1.0, 1.0}, Face_t::POSZ, 1.73205080756888, c},
            {{0.0, -1.0, -1.0}, Face_t::NEGZ, 1.41421356237310, {2, 3, 2}},
            {{0.0, -1.0, 1.0}, Face_t::POSZ, 1.41421356237310, c},
            {{0.0, 1.0, -1.0}, Face_t::NEGZ, 1.41421356237310, {2, 3, 2}},
            {{0.0, 1.0, 1.0}, Face_t::POSZ, 1.41421356237310, c},
            {{1.0, -1.0, -1.0}, Face_t::NEGZ, 1.73205080756888, {2, 3, 2}},
            {{1.0, -1.0, 0.0}, Face_t::NEG_W, 10.35276180410083, c},
            {{1.0, -1.0, 1.0}, Face_t::POSZ, 1.73205080756888, c},
            {{1.0, 0.0, -1.0}, Face_t::NEGZ, 1.41421356237310, {2, 3, 2}},
            {{1.0, 0.0, 0.0}, Face_t::POS_U, 10.00000000000000, c},
            {{1.0, 0.0, 1.0}, Face_t::POSZ, 1.41421356237310, c},
            {{1.0, 1.0, -1.0}, Face_t::NEGZ, 1.73205080756888, {2, 3, 2}},
            {{1.0, 1.0, 0.0}, Face_t::POS_V, 10.35276180410083, {2, 4, 3}},
            {{1.0, 1.0, 1.0}, Face_t::POSZ, 1.73205080756888, c},
            {{0.1, 0.1, 1.0}, Face_t::POSZ, 1.00995049383621, c},
            {{0.1, 0.1, -1.0}, Face_t::NEGZ, 1.00995049383621, {2, 3, 2}}};

        this->test_intersect(this->centroid(c), c, make_span(tests));
    }
}
