//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstConeAligned.cc
 * \brief ConeAligned class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../ConeAligned.hh"

#include "celeritas_test.hh"
#include "base/VectorFunctions.hh"

using namespace celeritas;
using celeritas::Real3;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(ConeXAligned, sense)
{
    ConeX cone(Real3{0, 0, 1}, 0.5);
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 1), cone.origin());
    EXPECT_SOFT_EQ(0.5 * 0.5, cone.tangent_sq());

    EXPECT_EQ(SignedSense::inside, cone.calc_sense(Real3{2, 0, 1}));
    EXPECT_EQ(SignedSense::inside, cone.calc_sense(Real3{2, 0.99, 1}));
    EXPECT_EQ(SignedSense::outside, cone.calc_sense(Real3{2, 1.01, 1}));

    EXPECT_EQ(SignedSense::inside, cone.calc_sense(Real3{2, 0, 1}));
    EXPECT_EQ(SignedSense::inside, cone.calc_sense(Real3{2, -0.99, 1}));
    EXPECT_EQ(SignedSense::outside, cone.calc_sense(Real3{2, -1.01, 1}));

    EXPECT_EQ(SignedSense::inside,
              cone.calc_sense(Real3{
                  -4,
                  0,
                  1,
              }));
    EXPECT_EQ(SignedSense::inside,
              cone.calc_sense(Real3{
                  -4,
                  0,
                  -0.99,
              }));
    EXPECT_EQ(SignedSense::outside,
              cone.calc_sense(Real3{
                  -4,
                  0,
                  -1.01,
              }));
}

TEST(ConeXAligned, normal)
{
    ConeX cone(Real3{2, 3, 4}, 0.5);

    Real3 pos(2 + 2, 3 + 1, 4 + 0);

    Real3 expected{-1, 2, 0};
    normalize_direction(&expected);
    EXPECT_VEC_SOFT_EQ(expected, cone.calc_normal(pos));
}

//---------------------------------------------------------------------------//
TEST(ConeXAligned, intersection_typical)
{
    // Cone with origin at (1.1, 2.2, 3.3) , slope of +- 2/3
    ConeX     cone(Real3(1.1, 2.2, 3.3), 2. / 3.);
    real_type distance[ConeX::num_intersections()] = {-1, -1};

    cone.calc_intersections(Real3(1.1 + 3, 2.2 + 2 + 1.0, 3.3),
                            Real3(0.0, -1.0, 0.0),
                            SurfaceState::off,
                            distance);
    EXPECT_SOFT_EQ(1.0, distance[0]);
    EXPECT_SOFT_EQ(5.0, distance[1]);
}

TEST(ConeXAligned, intersection_on_surface)
{
    // Cone with origin at (1.1, 2.2, 3.3) , slope of +- 2/3
    ConeX     cone(Real3(1.1, 2.2, 3.3), 2. / 3.);
    real_type distance[ConeX::num_intersections()] = {-1, -1};

    cone.calc_intersections(Real3(1.1 + 3, 2.2 + 2, 3.3),
                            Real3(0.0, -1.0, 0.0),
                            SurfaceState::on,
                            distance);
    EXPECT_SOFT_EQ(4.0, distance[0]);
    EXPECT_SOFT_EQ(no_intersection(), distance[1]);
}

TEST(ConeXAligned, intersection_along_surface)
{
    // Cone with origin at (1.1, 2.2, 3.3) , slope of +- 2/3
    ConeX     cone(Real3(1.1, 2.2, 3.3), 2. / 3.);
    real_type distance[ConeX::num_intersections()] = {-1, -1};

    // Along the cone edge heading up and right
    Real3 dir{3.0, 2.0, 0.0};
    normalize_direction(&dir);

    // below lower left sheet
    cone.calc_intersections(
        Real3(1.1 - 3, 2.2 - 2 - 1, 3.3), dir, SurfaceState::off, distance);
    EXPECT_SOFT_EQ(4.5069390943299865, distance[0]);
    EXPECT_SOFT_EQ(no_intersection(), distance[1]);

    // inside left sheet
    cone.calc_intersections(
        Real3(1.1 - 3, 2.2 - 2 + 1, 3.3), dir, SurfaceState::off, distance);
    EXPECT_SOFT_EQ(2.7041634565979917, distance[0]);
    EXPECT_SOFT_EQ(no_intersection(), distance[1]);

    // inside right sheet
    cone.calc_intersections(
        Real3(1.1 + 3, 2.2 + 2 - 1, 3.3), dir, SurfaceState::off, distance);
    EXPECT_SOFT_EQ(no_intersection(), distance[0]);
    EXPECT_SOFT_EQ(no_intersection(), distance[1]);

    // outside right sheet
    cone.calc_intersections(
        Real3(1.1 + 3, 2.2 + 2 + 1, 3.3), dir, SurfaceState::off, distance);
    EXPECT_SOFT_EQ(no_intersection(), distance[0]);
    EXPECT_SOFT_EQ(no_intersection(), distance[1]);

    // DOUBLE DEGENERATE: on sheet, traveling along it
    cone.calc_intersections(
        Real3(1.1 + 3, 2.2 + 2, 3.3), dir, SurfaceState::on, distance);
    EXPECT_SOFT_EQ(no_intersection(), distance[0]);
    EXPECT_SOFT_EQ(no_intersection(), distance[1]);
}

//---------------------------------------------------------------------------//
TEST(ConeXAligned, debug_intersect)
{
    ConeX cone(Real3(6.5, 4.0, 0.0), 2. / 3.);
    Real3 pos(11.209772802620732, 4.7633848351978276, 1.438286089005397);
    Real3 dir(
        -0.97544487122446721, -0.21994072746270818, -0.011549008834416619);

    real_type distance[ConeX::num_intersections()] = {-1, -1};
    cone.calc_intersections(pos, dir, SurfaceState::off, distance);
    EXPECT_SOFT_EQ(2.645662863301006, distance[0]);
    EXPECT_SOFT_EQ(9.9221689631064418 - 2.645662863301006, distance[1]);

    pos += dir * distance[0];
    real_type expected_next = distance[1] - distance[0];
    cone.calc_intersections(pos, dir, SurfaceState::on, distance);
    EXPECT_SOFT_EQ(expected_next, distance[0]);
    EXPECT_SOFT_EQ(no_intersection(), distance[1]);

    pos += dir * distance[0];
    cone.calc_intersections(pos, dir, SurfaceState::on, distance);
    EXPECT_SOFT_EQ(no_intersection(), distance[0]);
    EXPECT_SOFT_EQ(no_intersection(), distance[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Test initialization on or near boundary
 */
class DegenerateBoundaryTest : public Test
{
  protected:
    void SetUp()
    {
        origin = {1.1, 2.2, 3.3};
        radius = 0.9;
        eps    = 0.0;
        z      = 1.0;
    }

    void run() const;

  protected:
    Real3     origin;
    real_type z;
    real_type radius;
    real_type eps;
};

void DegenerateBoundaryTest::run() const
{
    ConeZ           cone(origin, radius);
    real_type       distances[] = {-1, -1, -1, -1};
    const real_type tol         = std::max(1.e-14, std::fabs(eps));

    // Distance across the cone at the current point
    const real_type diameter = 2 * std::fabs(z) * radius;

    Real3 pos = origin;
    // move so that the circular cross section looks like an equivalent
    // cylinder
    pos[2] += 1.0;
    Real3 dir = {1, 0, 0};

    //// Left boundary ////
    pos[0] = origin[0] - diameter / 2 - eps;

    cone.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_NEAR(diameter + eps, distances[0], tol);
    EXPECT_EQ(no_intersection(), distances[1]);

    //// Right boundary ////
    pos[0] = origin[0] + diameter / 2 + eps;

    cone.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // Make sure we never go off the end
    EXPECT_EQ(-1.0, distances[2]);
    EXPECT_EQ(-1.0, distances[3]);
}

TEST_F(DegenerateBoundaryTest, neg_above)
{
    z   = 1.0;
    eps = -1.0e-8;
    run();
}

TEST_F(DegenerateBoundaryTest, zero_above)
{
    z   = 1.0;
    eps = 0.0;
    run();
}

TEST_F(DegenerateBoundaryTest, pos_above)
{
    z   = 1.0;
    eps = 1.e-8;
    run();
}

TEST_F(DegenerateBoundaryTest, neg_below)
{
    z   = -1.0;
    eps = -1.0e-8;
    run();
}

TEST_F(DegenerateBoundaryTest, zero_below)
{
    z   = -1.0;
    eps = 0.0;
    run();
}

TEST_F(DegenerateBoundaryTest, pos_below)
{
    z   = -1.0;
    eps = 1.e-8;
    run();
}
