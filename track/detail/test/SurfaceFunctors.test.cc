//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/test/tstSurfaceFunctors.cc
 * \brief Tests for class SurfaceFunctors
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../SurfaceFunctors.hh"

#include "celeritas_test.hh"
#include "base/FixedViewArray.hh"
#include "orange/surfaces/PlaneAligned.hh"
#include "orange/surfaces/Sphere.hh"
#include "orange/surfaces/SurfaceAction.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class SurfaceFunctorsTest : public ::Test
{
  protected:
    //// TYPE ALIASES ////

  protected:
    void SetUp()
    {
        surfaces.push_back(PlaneX(1.25));
        surfaces.push_back(Sphere({2.25, 1, 0}, 1.25));
    }

  protected:
    //// DATA ////

    celeritas::SurfaceContainer surfaces;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(SurfaceFunctorsTest, calc_sense)
{
    using celeritas::detail::CalcSense;
    Real3     pos{0.9, 0, 0};
    CalcSense calc{make_fixed_view(pos)};

    EXPECT_EQ(SignedSense::inside, calc(surfaces.get<PlaneX>(SurfaceId{0})));
    EXPECT_EQ(SignedSense::outside, calc(surfaces.get<Sphere>(SurfaceId{1})));

    pos = {1.0, 0, 0};
    EXPECT_EQ(SignedSense::inside, calc(surfaces.get<PlaneX>(SurfaceId{0})));
    EXPECT_EQ(SignedSense::outside, calc(surfaces.get<Sphere>(SurfaceId{1})));

    // Test as generic surfaces
    pos = {2, 0, 0};
    auto calc_generic
        = make_surface_action(surfaces, CalcSense{make_fixed_view(pos)});
    EXPECT_EQ(SignedSense::outside, calc_generic(SurfaceId{0}));
    EXPECT_EQ(SignedSense::inside, calc_generic(SurfaceId{1}));
}

//---------------------------------------------------------------------------//

TEST_F(SurfaceFunctorsTest, num_intersections)
{
    auto num_intersections
        = make_surface_action(surfaces, celeritas::detail::NumIntersections{});
    EXPECT_EQ(1, num_intersections(SurfaceId{0}));
    EXPECT_EQ(2, num_intersections(SurfaceId{1}));
}

//---------------------------------------------------------------------------//

TEST_F(SurfaceFunctorsTest, calc_normal)
{
    Real3 pos;
    auto  calc_normal = make_surface_action(
        surfaces, celeritas::detail::CalcNormal{make_fixed_view(pos)});

    pos = {1.25, 1, 1};
    EXPECT_EQ(Real3({1, 0, 0}), calc_normal(SurfaceId{0}));
    pos = {2.25, 2.25, 0};
    EXPECT_EQ(Real3({0, 1, 0}), calc_normal(SurfaceId{1}));
}

//---------------------------------------------------------------------------//

TEST_F(SurfaceFunctorsTest, calc_safety_distance)
{
    Real3 pos;

    auto calc_distance = make_surface_action(
        surfaces, celeritas::detail::CalcSafetyDistance{make_fixed_view(pos)});

    real_type eps = 1e-4;
    pos           = {1.25 + eps, 1, 0};
    EXPECT_SOFT_EQ(eps, calc_distance(SurfaceId{0}));
    EXPECT_SOFT_EQ(0.25 + eps, calc_distance(SurfaceId{1}));

    pos = {1.25, 1, 0};
    EXPECT_SOFT_EQ(0, calc_distance(SurfaceId{0}));
    EXPECT_SOFT_EQ(0.25, calc_distance(SurfaceId{1}));

    pos = {1.25 - eps, 1, 0};
    EXPECT_SOFT_EQ(eps, calc_distance(SurfaceId{0}));
    EXPECT_SOFT_EQ(0.25 - eps, calc_distance(SurfaceId{1}));

    pos = {1.0 - eps, 1, 0};
    EXPECT_SOFT_EQ(0.25 + eps, calc_distance(SurfaceId{0}));
    EXPECT_SOFT_EQ(eps, calc_distance(SurfaceId{1}));

    pos = {3.5 + eps, 1, 0};
    EXPECT_SOFT_EQ(2.25 + eps, calc_distance(SurfaceId{0}));
    EXPECT_SOFT_NEAR(0.0 + eps, calc_distance(SurfaceId{1}), 1e-11);

    pos = {3.5, 1, 0};
    EXPECT_SOFT_EQ(2.25, calc_distance(SurfaceId{0}));
    EXPECT_SOFT_EQ(0.0, calc_distance(SurfaceId{1}));
}
