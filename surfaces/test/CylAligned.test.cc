//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstCylAligned.cc
 * \brief CylAligned class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../CylAligned.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/TransformUtils.hh"
#include "orange/Transform.hh"
#include "../GeneralQuadric.hh"
#include "../SimpleQuadric.hh"

using namespace geometria;
using namespace celeritas;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(CylXAligned, sense)
{
    CylX cyl(Real3{0, 0, 0}, 4.0);

    cout << cyl << "\n";

    EXPECT_EQ(SignedSense::inside, cyl.calc_sense(Real3{0, 3, 0}));
    EXPECT_EQ(SignedSense::outside, cyl.calc_sense(Real3{0, 0, 9}));
}

TEST(CylXAligned, normal)
{
    Real3 norm;

    CylX cyl(Real3{0, 1.23, 2.34}, 3.45);

    norm = cyl.calc_normal(Real3{1.23, 4.68, 2.34});
    EXPECT_SOFT_EQ(0, norm[0]);
    EXPECT_SOFT_EQ(1, norm[1]);
    EXPECT_SOFT_EQ(0, norm[2]);

    norm = cyl.calc_normal(Real3{12345., 2.7728869044748548, 5.42577380894971});
    EXPECT_SOFT_EQ(0, norm[0]);
    EXPECT_SOFT_EQ(0.4472135954999578, norm[1]);
    EXPECT_SOFT_EQ(0.894427190999916, norm[2]);
}

TEST(CylXAligned, intersect)
{
    real_type distances[] = {-1, -1};

    // From inside
    CylX cyl(Real3{1234.5, 0, 1}, 3.0);
    cyl.calc_intersections(
        Real3{0, 0, 2.5}, Real3{0, 1, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.598076211353316, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting both
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{0, 0, -4.0}, Real3{0, 0, 1}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.0, distances[0]);
    EXPECT_SOFT_EQ(8.0, distances[1]);

    // From outside, tangent
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{0, -3, -4.0}, Real3{0, 0, 1}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(5.0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting neither
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{0, 0, -4.0}, Real3{0, 1, 0}, SurfaceState::off, distances);

    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

TEST(CylXAligned, intersect_from_surface)
{
    real_type distances[CylX::num_intersections()];

    CylX cyl(Real3{0, 1.23, 2.34}, 3.45);

    // One intercept

    cyl.calc_intersections(
        Real3{1.23, 4.68, 2.34}, Real3{0, -1, 0}, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(6.9, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    cyl.calc_intersections(Real3{12345., 2.7728869044748548, 5.42577380894971},
                           Real3{0.6, 0, -0.8},
                           SurfaceState::on,
                           distances);
    EXPECT_SOFT_EQ(7.714434522374273, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // No intercepts

    cyl.calc_intersections(
        Real3{1.23, 4.68, 2.34}, Real3{0, 1, 0}, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    cyl.calc_intersections(Real3{12345., 2.7728869044748548, 5.42577380894971},
                           Real3{0.6, 0, 0.8},
                           SurfaceState::on,
                           distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

//---------------------------------------------------------------------------//

TEST(CylYAligned, sense)
{
    CylY cyl(Real3{1, 1234.5, 0}, 3.0);

    cout << cyl << "\n";

    EXPECT_EQ(SignedSense::inside, cyl.calc_sense(Real3{2.5, 0, 0}));
    EXPECT_EQ(SignedSense::outside, cyl.calc_sense(Real3{4.01, 0, 0}));
}

TEST(CylYAligned, normal)
{
    Real3 norm;

    CylY cyl(Real3{1.23, 0, 2.34}, 3.45);

    norm = cyl.calc_normal(Real3{4.68, 1.23, 2.34});
    EXPECT_SOFT_EQ(1, norm[0]);
    EXPECT_SOFT_EQ(0, norm[1]);
    EXPECT_SOFT_EQ(0, norm[2]);

    norm = cyl.calc_normal(Real3{2.7728869044748548, 12345., 5.42577380894971});
    EXPECT_SOFT_EQ(0.4472135954999578, norm[0]);
    EXPECT_SOFT_EQ(0, norm[1]);
    EXPECT_SOFT_EQ(0.894427190999916, norm[2]);
}

TEST(CylYAligned, intersect)
{
    real_type distances[] = {-1, -1};

    // From inside
    CylY cyl(Real3{1, 1234.5, 0}, 3.0);
    cyl.calc_intersections(
        Real3{2.5, 0, 0}, Real3{0, 0, 1}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.598076211353316, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting both
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-4.0, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.0, distances[0]);
    EXPECT_SOFT_EQ(8.0, distances[1]);

    // From outside, tangent
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-4.0, 0, -3}, Real3{1, 0, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(5.0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting neither
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-4.0, 0, 0}, Real3{0, 0, 1}, SurfaceState::off, distances);

    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

TEST(CylYAligned, intersect_from_surface)
{
    real_type distances[CylY::num_intersections()];

    CylY cyl(Real3{1.23, 0, 2.34}, 3.45);

    // One intercept

    cyl.calc_intersections(
        Real3{4.68, 1.23, 2.34}, Real3{-1, 0, 0}, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(6.9, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    cyl.calc_intersections(Real3{2.7728869044748548, 12345., 5.42577380894971},
                           Real3{0, 0.6, -0.8},
                           SurfaceState::on,
                           distances);
    EXPECT_SOFT_EQ(7.714434522374273, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // No intercepts

    cyl.calc_intersections(
        Real3{4.68, 1.23, 2.34}, Real3{1, 0, 0}, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    cyl.calc_intersections(Real3{2.7728869044748548, 12345., 5.42577380894971},
                           Real3{0, 0.6, 0.8},
                           SurfaceState::on,
                           distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

//---------------------------------------------------------------------------//

TEST(CylZAligned, sense)
{
    CylZ cyl(Real3{1, 0, 1234.5}, 3.0);

    cout << cyl << "\n";

    EXPECT_EQ(SignedSense::inside, cyl.calc_sense(Real3{2.5, 0, 0}));
    EXPECT_EQ(SignedSense::outside, cyl.calc_sense(Real3{4.01, 0, 0}));
}

TEST(CylZAligned, normal)
{
    Real3 norm;

    CylZ cyl(Real3{1.23, 2.34, 0}, 3.45);

    norm = cyl.calc_normal(Real3{4.68, 2.34, 1.23});
    EXPECT_SOFT_EQ(1, norm[0]);
    EXPECT_SOFT_EQ(0, norm[1]);
    EXPECT_SOFT_EQ(0, norm[2]);

    norm = cyl.calc_normal(Real3{2.7728869044748548, 5.42577380894971, 12345.});
    EXPECT_SOFT_EQ(0.4472135954999578, norm[0]);
    EXPECT_SOFT_EQ(0.894427190999916, norm[1]);
    EXPECT_SOFT_EQ(0, norm[2]);
}

TEST(CylZAligned, calc_intersections)
{
    real_type distances[] = {-1, -1};

    // From inside
    CylZ cyl(Real3{1, 0, 1234.5}, 3.0);
    cyl.calc_intersections(
        Real3{2.5, 0, 0}, Real3{0, 1, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.598076211353316, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting both
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-4.0, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.0, distances[0]);
    // TODO: we should calculate intersection to both surfaces
    // EXPECT_EQ(8.0,      distances[1]);

    // From outside, tangent
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-4.0, -3, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(5.0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting neither
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-4.0, 0, 0}, Real3{0, 1, 0}, SurfaceState::off, distances);

    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

TEST(CylZAligned, calc_intersections_on_surface)
{
    real_type distances[CylZ::num_intersections()];

    CylZ cyl(Real3{1.23, 2.34, 0}, 3.45);

    // One intercept

    cyl.calc_intersections(
        Real3{4.68, 2.34, 1.23}, Real3{-1, 0, 0}, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(6.9, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    cyl.calc_intersections(Real3{2.7728869044748548, 5.42577380894971, 12345.},
                           Real3{0, -0.8, 0.6},
                           SurfaceState::on,
                           distances);
    EXPECT_SOFT_EQ(7.714434522374273, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // No intercepts

    cyl.calc_intersections(
        Real3{4.68, 2.34, 1.23}, Real3{1, 0, 0}, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    cyl.calc_intersections(Real3{2.7728869044748548, 5.42577380894971, 12345.},
                           Real3{0, 0.8, 0.6},
                           SurfaceState::on,
                           distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

TEST(CylZAligned, calc_intersections_degenerate)
{
    real_type       distances[CylZ::num_intersections()];
    const real_type eps = 1.e-4;

    {
        CylZ cyl(Real3{1.0, 0.0, 0}, 1.0);

        // Heading toward, slightly inside
        cyl.calc_intersections(
            Real3{eps, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);
        EXPECT_SOFT_EQ(2.0 - eps, distances[0]);
        EXPECT_SOFT_EQ(no_intersection(), distances[1]);

        // Heading away, slightly inside
        cyl.calc_intersections(
            Real3{eps, 0, 0}, Real3{-1, 0, 0}, SurfaceState::off, distances);
        EXPECT_SOFT_EQ(eps, distances[0]);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Tangent, inside, off surface
        cyl.calc_intersections(
            Real3{eps, 0, 0}, Real3{0, 1, 0}, SurfaceState::off, distances);
        EXPECT_SOFT_NEAR(std::sqrt(2 * eps), distances[0], eps);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Tangent, inside, on surface
        cyl.calc_intersections(
            Real3{eps, 0, 0}, Real3{0, 1, 0}, SurfaceState::on, distances);
        EXPECT_EQ(no_intersection(), distances[0]);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Heading in
        cyl.calc_intersections(
            Real3{eps, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);
        EXPECT_SOFT_NEAR(2.0 + eps, distances[0], eps);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Heading away
        cyl.calc_intersections(Real3{2.0 - eps, 0, 0},
                               Real3{1, 0, 0},
                               SurfaceState::off,
                               distances);
        EXPECT_SOFT_NEAR(eps, distances[0], eps);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Tangent
        cyl.calc_intersections(
            Real3{eps, 0, 0}, Real3{0, 1, 0}, SurfaceState::off, distances);
        EXPECT_SOFT_EQ(std::sqrt(2 * eps * 1.0 - ipow<2>(eps)), distances[0]);
        EXPECT_EQ(no_intersection(), distances[1]);
    }
    {
        CylZ cyl(Real3{1.23, 2.34, 0}, 3.45);

        cyl.calc_intersections(Real3{4.68 - eps, 2.34, 1.23},
                               Real3{-1, 0, 0},
                               SurfaceState::off,
                               distances);
        EXPECT_SOFT_EQ(6.9 - eps, distances[0]);
        EXPECT_EQ(no_intersection(), distances[1]);

        cyl.calc_intersections(Real3{4.68 - eps, 2.34, 1.23},
                               Real3{1, 0, 0},
                               SurfaceState::off,
                               distances);
        EXPECT_SOFT_NEAR(eps, distances[0], eps);
        EXPECT_EQ(no_intersection(), distances[1]);
    }
}

TEST(CylZAligned, rotation)
{
    CylZ cyl(Real3{1, 2, 4}, 3.0);

    // Rotate 90 degrees clockwise about z
    {
        Transform t(rotation_matrix(Axis::z, .25));
        auto      gq = cyl.transformed(t);
        ASSERT_TRUE(SimpleQuadric::can_simplify(gq));
        SimpleQuadric sq(gq);
        ASSERT_TRUE(CylZ::can_simplify(sq));

        CylZ rc(sq);
        EXPECT_SOFT_EQ(-2, rc.origin_u());
        EXPECT_SOFT_EQ(1, rc.origin_v());
        EXPECT_SOFT_EQ(3.0 * 3.0, rc.radius_sq());
    }

    // Rotate 180 degrees about z
    {
        Transform t(rotation_matrix(Axis::z, .5));
        auto      gq = cyl.transformed(t);
        ASSERT_TRUE(SimpleQuadric::can_simplify(gq));
        SimpleQuadric sq(gq);
        ASSERT_TRUE(CylZ::can_simplify(sq));

        CylZ rc(sq);
        EXPECT_SOFT_EQ(-1, rc.origin_u());
        EXPECT_SOFT_EQ(-2, rc.origin_v());
        EXPECT_SOFT_EQ(3.0 * 3.0, rc.radius_sq());
    }

    // Rotate 90 degrees about x
    {
        Transform t(rotation_matrix(Axis::x, .25));
        auto      gq = cyl.transformed(t);
        ASSERT_TRUE(SimpleQuadric::can_simplify(gq));
        SimpleQuadric sq(gq);
        ASSERT_TRUE(CylY::can_simplify(sq));

        CylY rc(sq);
        EXPECT_SOFT_EQ(1, rc.origin_u());
        EXPECT_SOFT_EQ(2, rc.origin_v());
        EXPECT_SOFT_EQ(3.0 * 3.0, rc.radius_sq());
    }

    // Rotate 90 degrees about y
    {
        Transform t(rotation_matrix(Axis::y, .25));
        auto      gq = cyl.transformed(t);
        ASSERT_TRUE(SimpleQuadric::can_simplify(gq));
        SimpleQuadric sq(gq);
        ASSERT_TRUE(CylX::can_simplify(sq));

        CylX rc(sq);
        EXPECT_SOFT_EQ(2, rc.origin_u());
        EXPECT_SOFT_EQ(-1, rc.origin_v());
        EXPECT_SOFT_EQ(3.0 * 3.0, rc.radius_sq());
    }

    // Rotate 45 degrees about x
    {
        Transform t(rotation_matrix(Axis::x, .5 / 4));
        auto      gq = cyl.transformed(t);
        EXPECT_FALSE(SimpleQuadric::can_simplify(gq));

        const real_type sqrt8 = std::sqrt(8.);
        EXPECT_VEC_SOFT_EQ(Real3(1, .5, .5), gq.second());
        EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), gq.cross());
        EXPECT_VEC_SOFT_EQ(Real3(-2, -sqrt8, -sqrt8), gq.first());
        EXPECT_SOFT_EQ(-4, gq.zeroth());
    }
}

/*! Note about disabled tests on exact cylinder boundary.
 *
 * In these tests, the intersection distance (in the two-real-root case of
 * solve_quadratic in Utils.cc) evaluates to exactly zero, which
 * is not returned as a viable solution (since there are other cases like when
 * a particle is near a surface but outside and headed away).
 */

// This is basically the same test as UnityDegenerateBoundaryTest.zero
TEST(CylZAligned, DISABLED_coincident_inside)
{
    CylZ      cyl(Real3{0, 0, 0}, 1);
    real_type distances[CylZ::num_intersections()];

    cyl.calc_intersections(
        Real3{1, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

// This is basically the same test as UnityDegenerateBoundaryTest.zero
TEST(CylZAligned, DISABLED_coincident_inside2)
{
    CylZ      cyl(Real3{0, 0, 0}, 0.5);
    real_type distances[CylZ::num_intersections()];

    cyl.calc_intersections(Real3{0.49999999932758, -2.5931062327113027e-05, 0},
                           Real3{1, 0, 0},
                           SurfaceState::off,
                           distances);
    EXPECT_SOFT_EQ(0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Test initialization on or near boundary
 *
 * These degenerate tests are no longer very useful since we do not require an
 * "inside" or "outside" sense to be passed to the cylinder.
 */
class DegenerateBoundaryTest : public Test
{
  protected:
    virtual void SetUp()
    {
        origin = {1.1, 2.2, 3.3};
        radius = 0.9;
        eps    = 0.0;
    }

    void run(real_type xdir) const;

  protected:
    Real3     origin;
    real_type radius;
    real_type eps;
};

void DegenerateBoundaryTest::run(real_type xdir) const
{
    CylZ            cyl(origin, radius);
    real_type       distances[] = {-1, -1, -1, -1};
    const real_type tol         = std::max(1.e-14, 2 * std::fabs(eps));

    const char* dirstr = (xdir == 1 ? "toward the right" : "toward the left");
    // Distance across the cylinder
    const real_type diameter = 2 * radius;

    Real3 pos = origin;
    Real3 dir = {xdir, 0, 0};
    SCOPED_TRACE(dirstr);

    //// Inward boundary ////
    pos[0] = origin[0] - xdir * (diameter / 2 + eps);
    cyl.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_NEAR(diameter + eps, distances[0], tol);
    EXPECT_EQ(no_intersection(), distances[1]);

    //// Outward boundary ////
    pos[0] = origin[0] + xdir * (diameter / 2 + eps);
    cyl.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // Make sure we never go off the end
    EXPECT_EQ(-1.0, distances[2]);
    EXPECT_EQ(-1.0, distances[3]);
}

TEST_F(DegenerateBoundaryTest, neg)
{
    eps = -1.0e-8;
    run(1);
    run(-1);
}

TEST_F(DegenerateBoundaryTest, zero)
{
    eps = 0.0;
    run(1);
    run(-1);
}

TEST_F(DegenerateBoundaryTest, pos)
{
    eps = 1.e-8;
    run(1);
    run(-1);
}

//---------------------------------------------------------------------------//
/*!
 * Test initialization on or near boundary
 */
class UnityDegenerateBoundaryTest : public DegenerateBoundaryTest
{
  protected:
    virtual void SetUp()
    {
        origin = {0, 0, 0};
        radius = 1.0;
        eps    = 0;
    }
};

TEST_F(UnityDegenerateBoundaryTest, neg)
{
    eps = -1.0e-8;
    run(1);
    run(-1);
}

TEST_F(UnityDegenerateBoundaryTest, DISABLED_zero)
{
    eps = 0.0;
    run(1);
    run(-1);
}

TEST_F(UnityDegenerateBoundaryTest, pos)
{
    eps = 1.e-8;
    run(1);
    run(-1);
}
