//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstPlaneAligned.cc
 * \brief PlaneAligned class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../PlaneAligned.hh"

#include <cmath>
#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/TransformUtils.hh"
#include "orange/Transform.hh"
#include "../Plane.hh"

using celeritas::no_intersection;
using celeritas::PlaneX;
using celeritas::PlaneY;
using celeritas::PlaneZ;
using celeritas::Real3;
using celeritas::SignedSense;
using celeritas::SurfaceState;

//---------------------------------------------------------------------------//
// FIXTURE
class PlaneTestAligned : public Test
{
  protected:
    template<class S>
    real_type
    calc_intersection(const S& surf, Real3 pos, Real3 dir, SurfaceState s)
    {
        real_type distances[] = {-1, -1, -1, -1};
        surf.calc_intersections(pos, dir, s, distances);

        // Make sure the surface hasn't modified the higher distances
        EXPECT_EQ(-1., distances[1]);
        EXPECT_EQ(-1., distances[2]);
        EXPECT_EQ(-1., distances[3]);
        return distances[0];
    }

    template<class S>
    real_type
    calc_intersection_from_surface(const S& surf, Real3 pos, Real3 dir)
    {
        real_type distances[] = {-1, -1, -1};
        surf.calc_intersections(pos, dir, SurfaceState::on, distances);

        EXPECT_EQ(-1., distances[1]);
        EXPECT_EQ(-1., distances[2]);
        return distances[0];
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(PlaneTestAligned, x_plane)
{
    PlaneX p(1.0);

    EXPECT_EQ(Real3(1.0, 0.0, 0.0), p.calc_normal(Real3(1.0, 0.1, -4.2)));

    // Test sense
    EXPECT_EQ(SignedSense::outside, p.calc_sense(Real3({1.01, 0.1, 0.0})));
    EXPECT_EQ(SignedSense::inside, p.calc_sense(Real3({0.99, 0.1, 0.0})));

    // Simple intersections
    {
        // Direction vectors
        Real3 px({1.0, 0.0, 0.0});
        Real3 mx({-1.0, 0.0, 0.0});
        Real3 py({0.0, 1.0, 0.0});
        Real3 pz({0.0, 0.0, 1.0});

        // On surface
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {0.9999, 0.0, 0.0}, px, SurfaceState::on));
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {1.0001, 0.0, 0.0}, mx, SurfaceState::on));

        // Intersections
        EXPECT_SOFT_EQ(
            0.0, calc_intersection(p, {1.0, 0.0, 0.0}, px, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.0, calc_intersection(p, {1.0, 0.0, 0.0}, mx, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.01,
            calc_intersection(p, {0.99, 0.0, 0.0}, px, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.01,
            calc_intersection(p, {1.01, 0.0, 0.0}, mx, SurfaceState::off));

        // parallel
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {1.01, 0.0, 0.0}, py, SurfaceState::off));
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {0.99, 0.0, 0.0}, pz, SurfaceState::off));
    }
    // Point on surface
    {
        real_type dist;

        Real3 xpos{1.0, 0.0, 0.0};

        dist = calc_intersection_from_surface(p, {1.0, 0.0, 0.0}, xpos);
        EXPECT_EQ(no_intersection(), dist);
    }
}

//---------------------------------------------------------------------------//

TEST_F(PlaneTestAligned, y_plane)
{
    PlaneY p(-1.0);

    EXPECT_EQ(Real3({0.0, 1.0, 0.0}), p.calc_normal(Real3{1.0, -1.0, -4.2}));

    // Test sense
    EXPECT_EQ(SignedSense::outside, p.calc_sense(Real3({1.01, -0.99, 0.0})));
    EXPECT_EQ(SignedSense::inside, p.calc_sense(Real3({0.99, -1.01, 0.0})));

    // Simple intersections
    {
        Real3 py({0.0, 1.0, 0.0});
        Real3 my({0.0, -1.0, 0.0});
        Real3 px({1.0, 0.0, 0.0});
        Real3 pz({0.0, 0.0, 1.0});

        // positive sense
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {0.0, -1.01, 0.0}, py, SurfaceState::on));

        // negative sense
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {0.0, -0.99, 0.0}, my, SurfaceState::on));

        // Intersections
        EXPECT_SOFT_EQ(
            0.0, calc_intersection(p, {0.0, -1.0, 0.0}, py, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.0, calc_intersection(p, {0.0, -1.0, 0.0}, my, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.01,
            calc_intersection(p, {0.0, -1.01, 0.0}, py, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.01,
            calc_intersection(p, {0.0, -0.99, 0.0}, my, SurfaceState::off));

        // parallel
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {-1.01, 1.0, 0.0}, px, SurfaceState::off));
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {-0.99, -1.1, 0.0}, pz, SurfaceState::off));
    }
    // Point on surface
    {
        real_type dist;

        Real3 yneg{0.0, -1.0, 0.0};

        dist = calc_intersection_from_surface(p, {0.0, 1.0, 0.0}, yneg);
        EXPECT_EQ(no_intersection(), dist);
    }
}

//---------------------------------------------------------------------------//

TEST_F(PlaneTestAligned, plane_z)
{
    PlaneZ p(0.0);

    EXPECT_EQ(Real3({0.0, 0.0, 1.0}), p.calc_normal(Real3{1.0, 0.1, 0.0}));

    // Test sense
    EXPECT_EQ(SignedSense::outside, p.calc_sense(Real3({1.01, 0.1, 0.01})));
    EXPECT_EQ(SignedSense::inside, p.calc_sense(Real3({0.99, 0.1, -0.01})));

    // Simple intersections
    {
        Real3 pz({0.0, 0.0, 1.0});
        Real3 mz({0.0, 0.0, -1.0});
        Real3 px({1.0, 0.0, 0.0});
        Real3 py({0.0, 1.0, 0.0});

        // positive sense
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {0.0, 0.0, -0.01}, pz, SurfaceState::on));

        // negative sense
        EXPECT_EQ(no_intersection(),
                  calc_intersection(p, {0.0, 0.0, 0.01}, mz, SurfaceState::on));

        // Intersections
        EXPECT_SOFT_EQ(
            0.0, calc_intersection(p, {0.0, 0.0, 0.0}, pz, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.0, calc_intersection(p, {0.0, 0.0, 0.0}, mz, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.01,
            calc_intersection(p, {0.0, 0.0, -0.01}, pz, SurfaceState::off));
        EXPECT_SOFT_EQ(
            0.01,
            calc_intersection(p, {0.0, 0.0, 0.01}, mz, SurfaceState::off));

        // parallel
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {-1.01, 0.0, 0.1}, px, SurfaceState::off));
        EXPECT_EQ(
            no_intersection(),
            calc_intersection(p, {-0.99, 0.0, -0.0}, py, SurfaceState::off));
    }
    // Point on surface
    {
        real_type dist;

        Real3 zpos{0.0, 0.0, 1.0};

        dist = calc_intersection_from_surface(p, {0.0, 0.0, 1.0}, zpos);
        EXPECT_EQ(no_intersection(), dist);
    }
}

//---------------------------------------------------------------------------//

TEST_F(PlaneTestAligned, rotation)
{
    using Axis::z;
    using celeritas::Transform;
    using geometria::rotation_matrix;

    PlaneX p(1.25);
    EXPECT_SOFT_EQ(1.25, p.position());

    // Rotate a quarter turn clockwise
    {
        auto rp = p.transformed(rotation_matrix(Z, .25));
        EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), rp.normal());
        EXPECT_SOFT_EQ(1.25, rp.displacement());

        EXPECT_FALSE(PlaneX::can_simplify(rp));
        EXPECT_TRUE(PlaneY::can_simplify(rp));
        EXPECT_FALSE(PlaneZ::can_simplify(rp));

        // Convert to y plane
        PlaneY yp(rp);
        EXPECT_SOFT_EQ(1.25, yp.position());
    }

    // Rotate a teeny bit more than a quarter turn clockwise: should still be
    // able to simplify
    {
        auto rp = p.transformed(rotation_matrix(Z, .25 + 1.e-11));

        EXPECT_FALSE(PlaneX::can_simplify(rp));
        EXPECT_TRUE(PlaneY::can_simplify(rp));
        EXPECT_FALSE(PlaneZ::can_simplify(rp));

        // Convert to y plane
        PlaneY yp(rp);
        EXPECT_SOFT_EQ(1.25, yp.position());
    }

    // Rotate so its normal is backward (half turn)
    {
        auto rp = p.transformed(rotation_matrix(Z, .5));
        EXPECT_VEC_SOFT_EQ(Real3(-1, 0, 0), rp.normal());
        EXPECT_SOFT_EQ(1.25, rp.displacement());

        // Can't simplify because the sense would change
        EXPECT_FALSE(PlaneX::can_simplify(rp));
        EXPECT_FALSE(PlaneY::can_simplify(rp));
        EXPECT_FALSE(PlaneZ::can_simplify(rp));
    }
}
