//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstSimpleQuadric.cc
 * \brief SimpleQuadric class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../SimpleQuadric.hh"

#include "celeritas_test.hh"
#include "orange/Transform.hh"
#include "orange/TransformUtils.hh"
#include "../GeneralQuadric.hh"

using namespace geometria;
using namespace celeritas;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Hyperboloid, print)
{
    /* From the MCNP manual example:
     *
     * This describes a surface symmetric about the x–axis, which passes
     * through the three (x,r) points (7,5), (3,2), and (4,3). This surface is
     * a hyperboloid of two sheets
     */
    SimpleQuadric sq(
        {-0.083333333, 1, 1}, {0, 0, 0}, 68.52083, Real3(-26.5, 0, 0));

    cout << sq << "\n";
}

TEST(Ellipsoid, origin)
{
    // 1 x 2.5 x .3 radii
    const Real3 second(
        2.5 * 2.5 * 0.3 * 0.3, 1.0 * 1.0 * 0.3 * 0.3, 1.0 * 1.0 * 2.5 * 2.5);
    const Real3     first(0, 0, 0);
    const real_type zeroth = -1 * 2.5 * 2.5 * 0.3 * 0.3;
    SimpleQuadric   sq(second, first, zeroth, {0, 0, 0});

    EXPECT_VEC_SOFT_EQ(second, sq.second());
    EXPECT_VEC_SOFT_EQ(first, sq.first());
    EXPECT_SOFT_EQ(zeroth, sq.zeroth());

    // Test intersections along major axes
    real_type dist[2] = {-1, -1};
    sq.calc_intersections(
        Real3(-2.5, 0, 0), Real3(1, 0, 0), SurfaceState::off, dist);
    EXPECT_SOFT_EQ(1.5, dist[0]);
    EXPECT_SOFT_EQ(1.5 + 2.0, dist[1]);
    sq.calc_intersections(
        Real3(0, 2.5, 0), Real3(0, -1, 0), SurfaceState::on, dist);
    EXPECT_SOFT_EQ(5.0, dist[0]);
    EXPECT_SOFT_EQ(no_intersection(), dist[1]);
    sq.calc_intersections(
        Real3(0, 0, 0), Real3(0, 0, 1), SurfaceState::off, dist);
    EXPECT_SOFT_EQ(0.3, dist[0]);
    EXPECT_SOFT_EQ(no_intersection(), dist[1]);
}

//---------------------------------------------------------------------------//

TEST(Ellipsoid, translation)
{
    SimpleQuadric orig({0.5625, .09, 6.25}, {0, 0, 0}, -0.5625, {0, 0, 0});

    SimpleQuadric sq = orig.translated(Transform(Real3(0.5, 4, 8)));

    // Test that it was constructed correctly
    {
        SimpleQuadric equiv(
            {0.5625, .09, 6.25}, {0, 0, 0}, -0.5625, {0.5, 4, 8});
        EXPECT_VEC_SOFT_EQ(equiv.second(), sq.second());
        EXPECT_VEC_SOFT_EQ(equiv.first(), sq.first());
        EXPECT_SOFT_EQ(equiv.zeroth(), sq.zeroth());
    }

    // Check quadric
    EXPECT_VEC_SOFT_EQ(Real3(0.5625, 0.09, 6.25), sq.second());
    EXPECT_VEC_SOFT_EQ(Real3(-0.5625, -0.72, -100), sq.first());
    EXPECT_SOFT_EQ(401.018125, sq.zeroth());

    // Test intersections along major axes
    real_type dist[2] = {-1, -1};
    sq.calc_intersections(
        Real3(-2.0, 4, 8), Real3(1, 0, 0), SurfaceState::off, dist);
    EXPECT_SOFT_EQ(1.5, dist[0]);
    EXPECT_SOFT_EQ(1.5 + 2.0, dist[1]);
    sq.calc_intersections(
        Real3(0.5, 7.5, 8.0), Real3(0, -1, 0), SurfaceState::off, dist);
    EXPECT_SOFT_EQ(1., dist[0]);
    EXPECT_SOFT_EQ(1 + 5.0, dist[1]);
    sq.calc_intersections(
        Real3(0.5, 4., 8.0), Real3(0, 0, -1), SurfaceState::off, dist);
    EXPECT_SOFT_EQ(0.3, dist[0]);
    EXPECT_SOFT_EQ(no_intersection(), dist[1]);
}

//---------------------------------------------------------------------------//

TEST(Ellipsoid, rotation)
{
    SimpleQuadric orig({0.5625, .09, 6.25}, {0, 0, 0}, -0.5625, {0, 0, 0});

    // Translate, then rotate a quarter turn about Z: x -> y, y->-x
    Transform t(rotation_matrix(Axis::z, .25), Real3(0.5, 4, 8));
    auto      gq = orig.transformed(t);

    ASSERT_TRUE(SimpleQuadric::can_simplify(gq));

    // Construct from rotated general quadric
    const SimpleQuadric sq(gq);

    // Check quadric
    EXPECT_VEC_SOFT_EQ(Real3(0.09, 0.5625, 6.25), sq.second());
    EXPECT_VEC_SOFT_EQ(Real3(-0.09, -4.5, -100), sq.first());
    EXPECT_SOFT_EQ(408.46, sq.zeroth());
}
