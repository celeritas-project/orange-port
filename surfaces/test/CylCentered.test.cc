//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstCylCentered.cc
 * \brief Tests for class CylCentered
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../CylCentered.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "base/ViewField.hh"
#include "orange/TransformUtils.hh"
#include "orange/Transform.hh"
#include "../GeneralQuadric.hh"
#include "../SimpleQuadric.hh"
#include "../CylAligned.hh"

using namespace geometria;
using namespace celeritas;
using cVFDbl = ViewField<const real_type>;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(CentXAligned, construction)
{
    EXPECT_EQ(1, CCylX::size());
    EXPECT_EQ(2, CCylX::num_intersections());

    CCylX c(4.0);

    const real_type expected_data[] = {ipow<2>(4)};

    EXPECT_VEC_SOFT_EQ(expected_data, cVFDbl(c.data(), c.data() + c.size()));
}

TEST(CentXAligned, sense)
{
    CCylX cyl(4.0);

    cout << cyl << "\n";

    EXPECT_EQ(SignedSense::inside, cyl.calc_sense(Real3{0, 3, 0}));
    EXPECT_EQ(SignedSense::outside, cyl.calc_sense(Real3{0, 0, 9}));
}

TEST(CentXAligned, normal)
{
    CCylX cyl(3.45);

    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), cyl.calc_normal(Real3{1.23, 3.45, 0}));
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 1), cyl.calc_normal(Real3{0.0, 0.0, 4.45}));

    // Normal from center of cylinder is ill-defined but shouldn't raise an
    // error
    auto norm = cyl.calc_normal(Real3{10.0, 0.0, 0.0});
    EXPECT_TRUE(std::isnan(norm[0]));
    EXPECT_TRUE(std::isnan(norm[1]));
    EXPECT_TRUE(std::isnan(norm[2]));
}

TEST(CentXAligned, intersect)
{
    real_type distances[] = {-1, -1};

    // From inside
    CCylX cyl(3.0);
    cyl.calc_intersections(
        Real3{0, 0, 1.5}, Real3{0, 1, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.598076211353316, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting both
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{0, 0, -5.0}, Real3{0, 0, 1}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.0, distances[0]);
    EXPECT_SOFT_EQ(8.0, distances[1]);

    // From outside, tangent
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{0, -3, -5.0}, Real3{0, 0, 1}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(5.0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting neither
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{0, 0, -5.0}, Real3{0, 1, 0}, SurfaceState::off, distances);

    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

TEST(CentXAligned, intersect_from_surface)
{
    real_type distances[CCylX::num_intersections()];

    CCylX cyl(3.45);

    // One intercept

    cyl.calc_intersections(
        Real3{1.23, 3.45, 0.0}, Real3{0, -1, 0}, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(6.9, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // No intercepts
    cyl.calc_intersections(
        Real3{1.23, 4.68, 2.34}, Real3{0, 1, 0}, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

//---------------------------------------------------------------------------//

TEST(CentYAligned, sense)
{
    CCylY cyl(3.0);

    cout << cyl << "\n";

    EXPECT_EQ(SignedSense::inside, cyl.calc_sense(Real3{1.5, 0, 0}));
    EXPECT_EQ(SignedSense::outside, cyl.calc_sense(Real3{3.01, 0, 0}));
}

TEST(CentYAligned, intersect)
{
    real_type distances[] = {-1, -1};

    // From inside
    CCylY cyl(3.0);
    cyl.calc_intersections(
        Real3{1.5, 0, 0}, Real3{0, 0, 1}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.598076211353316, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting both
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-5.0, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.0, distances[0]);
    EXPECT_SOFT_EQ(8.0, distances[1]);

    // From outside, tangent
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-5.0, 0, -3}, Real3{1, 0, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(5.0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting neither
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-5.0, 0, 0}, Real3{0, 0, 1}, SurfaceState::off, distances);

    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

TEST(CentYAligned, intersect_from_surface)
{
    real_type distances[CCylY::num_intersections()];

    CCylY cyl(3.45);

    // One intercept

    cyl.calc_intersections(
        Real3{3.45, 1.23, 0}, Real3{-1, 0, 0}, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(6.9, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    cyl.calc_intersections(Real3{1.5528869044748548, 12345., 3.08577380894971},
                           Real3{0, 0.6, -0.8},
                           SurfaceState::on,
                           distances);
    EXPECT_SOFT_EQ(7.714434522374273, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // No intercepts

    cyl.calc_intersections(
        Real3{3.45, 1.23, 0.0}, Real3{1, 0, 0}, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    cyl.calc_intersections(Real3{1.5528869044748548, 12345., 3.08577380894971},
                           Real3{0, 0.6, 0.8},
                           SurfaceState::on,
                           distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

//---------------------------------------------------------------------------//

TEST(CentZAligned, sense)
{
    CCylZ cyl(3.0);

    cout << cyl << "\n";

    EXPECT_EQ(SignedSense::inside, cyl.calc_sense(Real3{1.5, 0, 0}));
    EXPECT_EQ(SignedSense::outside, cyl.calc_sense(Real3{3.01, 0, 0}));
}

TEST(CentZAligned, calc_intersections)
{
    real_type distances[] = {-1, -1};

    // From inside
    CCylZ cyl(3.0);
    cyl.calc_intersections(
        Real3{1.5, 0, 0}, Real3{0, 1, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.598076211353316, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting both
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-5.0, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(2.0, distances[0]);
    EXPECT_EQ(8.0, distances[1]);

    // From outside, tangent
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-5.0, -3, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);

    EXPECT_SOFT_EQ(5.0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // From outside, hitting neither
    distances[0] = distances[1] = -1;
    cyl.calc_intersections(
        Real3{-5.0, 0, 0}, Real3{0, 1, 0}, SurfaceState::off, distances);

    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

TEST(CentZAligned, calc_intersections_on_surface)
{
    real_type       distances[CCylZ::num_intersections()];
    const real_type eps = 1.e-4;

    {
        CCylZ cyl(1.0);

        // Heading toward, slightly inside
        cyl.calc_intersections(
            Real3{-1 + eps, 0, 0}, Real3{1, 0, 0}, SurfaceState::on, distances);
        EXPECT_SOFT_NEAR(2.0 - eps, distances[0], eps);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Heading away, slightly inside
        cyl.calc_intersections(
            Real3{-1 + eps, 0, 0}, Real3{-1, 0, 0}, SurfaceState::on, distances);
        EXPECT_EQ(no_intersection(), distances[0]);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Tangent
        cyl.calc_intersections(
            Real3{-1 + eps, 0, 0}, Real3{0, 1, 0}, SurfaceState::on, distances);
        EXPECT_EQ(no_intersection(), distances[0]);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Heading in, slightly inside
        cyl.calc_intersections(
            Real3{-1 + eps, 0, 0}, Real3{1, 0, 0}, SurfaceState::on, distances);
        EXPECT_SOFT_NEAR(2.0 - eps, distances[0], eps);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Heading away, slightly outside
        cyl.calc_intersections(
            Real3{1.0 - eps, 0, 0}, Real3{1, 0, 0}, SurfaceState::on, distances);
        EXPECT_EQ(no_intersection(), distances[0]);
        EXPECT_EQ(no_intersection(), distances[1]);

        // Tangent, slightly inside
        cyl.calc_intersections(
            Real3{-1 + eps, 0, 0}, Real3{0, 1, 0}, SurfaceState::on, distances);
        EXPECT_EQ(no_intersection(), distances[0]);
        EXPECT_EQ(no_intersection(), distances[1]);
    }
}

// Fused multiply-add on some CPUs in opt mode can cause the results of
// nearly-tangent cylinder checking to change.
TEST(CentZAligned, multi_intersect)
{
    using Axis::y;

    CCylZ cyl(10.0);

    std::vector<double> all_first_distances;
    std::vector<double> all_distances;
    std::vector<double> all_y;
    for (real_type x : {-10 + 1e-7, -1e-7, -0.0, 0.0, 1e-7, 10 - 1e-7})
    {
        for (real_type v : {-1.0, 1.0})
        {
            Real3 pos(x, -10.0001 * v, 0);
            Real3 dir(0, v, 0);

            real_type d[2];

            // Transport to inside of cylinder
            cyl.calc_intersections(pos, dir, SurfaceState::off, d);
            ASSERT_NE(no_intersection(), d[0]);
            all_first_distances.push_back(d[0]);
            pos[Y] += d[0] * dir[Y];
            all_y.push_back(pos[Y]);

            // Transport to other side of cylinder
            cyl.calc_intersections(pos, dir, SurfaceState::on, d);
            all_distances.push_back(d[0]);
            ASSERT_NE(no_intersection(), d[0]);
            pos[Y] += d[0] * dir[Y];

            // We're done
            cyl.calc_intersections(pos, dir, SurfaceState::on, d);
            EXPECT_EQ(no_intersection(), d[0]);
        }
    }

    // PRINT_EXPECTED(all_first_distances);
    const real_type expected_all_first_distances[] = {9.99869,
                                                      9.99869,
                                                      0.0001,
                                                      0.0001,
                                                      0.0001,
                                                      0.0001,
                                                      0.0001,
                                                      0.0001,
                                                      0.0001,
                                                      0.0001,
                                                      9.99869,
                                                      9.99869};
    EXPECT_VEC_NEAR(expected_all_first_distances, all_first_distances, 1e-5);

    // PRINT_EXPECTED(all_y);
    const real_type expected_all_y[] = {0.00141421,
                                        -0.00141421,
                                        10,
                                        -10,
                                        10,
                                        -10,
                                        10,
                                        -10,
                                        10,
                                        -10,
                                        0.00141421,
                                        -0.00141421};
    EXPECT_VEC_NEAR(expected_all_y, all_y, 1e-5);

    // PRINT_EXPECTED(all_distances);
    const real_type expected_all_distances[] = {0.00282843,
                                                0.00282843,
                                                20,
                                                20,
                                                20,
                                                20,
                                                20,
                                                20,
                                                20,
                                                20,
                                                0.00282843,
                                                0.00282843};
    EXPECT_VEC_NEAR(expected_all_distances, all_distances, 1e-5);
}

TEST(CentZAligned, translation)
{
    CCylZ cyl(3.0);

    Transform t(Real3(3, 2, 1));
    CylZ      oc = cyl.translated(t);
    EXPECT_FALSE(CCylZ::can_simplify(oc));

    EXPECT_SOFT_EQ(3, oc.origin_u());
    EXPECT_SOFT_EQ(2, oc.origin_v());
    EXPECT_SOFT_EQ(3.0 * 3.0, oc.radius_sq());
}

TEST(CentZAligned, rotation)
{
    CCylZ cyl(3.0);

    // Rotate 90 degrees clockwise about z
    {
        Transform t(rotation_matrix(Axis::z, .25));
        auto      gq = cyl.transformed(t);

        ASSERT_TRUE(SimpleQuadric::can_simplify(gq));
        SimpleQuadric sq(gq);

        ASSERT_TRUE(CylZ::can_simplify(sq));
        CylZ rc(sq);
        EXPECT_SOFT_EQ(0, rc.origin_u());
        EXPECT_SOFT_EQ(0, rc.origin_v());
        EXPECT_SOFT_EQ(3.0 * 3.0, rc.radius_sq());

        ASSERT_TRUE(CCylZ::can_simplify(rc));
        CCylZ coc(rc);
        EXPECT_SOFT_EQ(3.0 * 3.0, coc.radius_sq());
    }

    // Rotate 90 degrees about x
    {
        Transform t(rotation_matrix(Axis::x, .25));
        auto      gq = cyl.transformed(t);

        ASSERT_TRUE(SimpleQuadric::can_simplify(gq));
        SimpleQuadric sq(gq);

        ASSERT_TRUE(CylY::can_simplify(sq));
        CylY rc(sq);
        EXPECT_SOFT_EQ(0, rc.origin_u());
        EXPECT_SOFT_EQ(0, rc.origin_v());
        EXPECT_SOFT_EQ(3.0 * 3.0, rc.radius_sq());

        EXPECT_TRUE(CCylY::can_simplify(rc));
    }
}

//---------------------------------------------------------------------------//
/*!
 * Test initialization on or near boundary
 */
class DegenerateBoundaryTest : public Test
{
  protected:
    void run(real_type xdir) const;

    void run_all()
    {
        for (real_type r : {0.93, 1.0})
        {
            radius = r;
            for (real_type dir : {1, -1})
            {
                std::ostringstream msg;
                msg << "r=" << radius << ", dir_x=" << dir;
                SCOPED_TRACE(msg.str());

                run(1);
                run(-1);
            }
        }
    }

  protected:
    real_type radius = -1;
    real_type eps    = -1;
};

void DegenerateBoundaryTest::run(real_type xdir) const
{
    CCylZ           cyl(radius);
    real_type       distances[] = {-1, -1, -1, -1};
    const real_type tol         = std::max(1.e-14, 2 * std::fabs(eps));

    // Distance across the cylinder
    const real_type diameter = 2 * radius;

    Real3 pos = {0, 0, 0};
    Real3 dir = {xdir, 0, 0};

    //// Inward boundary ////
    pos[0] = -xdir * (diameter / 2 + eps);
    cyl.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_NEAR(diameter + eps, distances[0], tol);
    EXPECT_EQ(no_intersection(), distances[1]);

    //// Outward boundary ////
    pos[0] = xdir * (diameter / 2 + eps);
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
    run_all();
}

TEST_F(DegenerateBoundaryTest, DISABLED_zero)
{
    eps = 0.0;
    run_all();
}

TEST_F(DegenerateBoundaryTest, pos)
{
    eps = 1.e-8;
    run_all();
}
