//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstPlane.cc
 * \brief Plane class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Plane.hh"

#include <cmath>
#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "base/VectorFunctions.hh"
#include "orange/TransformUtils.hh"
#include "orange/Transform.hh"

using celeritas::no_intersection;
using celeritas::SignedSense;
using celeritas::SurfaceState;

//---------------------------------------------------------------------------//
// Test fixture
class PlaneTest : public Test
{
  protected:
    using Plane = celeritas::Plane;
    using Real3 = celeritas::Real3;

  protected:
    real_type
    calc_intersection(const Plane& surf, Real3 pos, Real3 dir, SurfaceState s)
    {
        real_type distances[] = {-1, -1, -1, -1};
        surf.calc_intersections(pos, dir, s, distances);

        // Make sure the surface hasn't modified the higher distances
        EXPECT_EQ(-1., distances[1]);
        EXPECT_EQ(-1., distances[2]);
        EXPECT_EQ(-1., distances[3]);
        return distances[0];
    }

    real_type
    calc_intersection_from_surface(const Plane& surf, Real3 pos, Real3 dir)
    {
        real_type distances[] = {-1, -1, -1};
        surf.calc_intersections(pos, dir, SurfaceState::on, distances);

        EXPECT_EQ(-1., distances[1]);
        EXPECT_EQ(-1., distances[2]);
        return distances[0];
    }
};

const real_type sqrt_two = std::sqrt(2);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(PlaneTest, user_construction)
{
    Real3     n(1.0, 2.0, 3.0);
    real_type d = 4.0 / vector_magnitude(n);
    normalize_direction(&n);
    Plane p(n, d);

    // Check the data
    const real_type* data = p.data();

    EXPECT_SOFT_NEAR(0.26726124, data[0], 1.e-7);
    EXPECT_SOFT_NEAR(0.53452248, data[1], 1.e-7);
    EXPECT_SOFT_NEAR(0.80178373, data[2], 1.e-7);
    EXPECT_SOFT_NEAR(1.06904497, data[3], 1.e-7);

    // Check accessors
    EXPECT_VEC_NEAR(
        Real3(0.26726124, 0.53452248, 0.80178373), p.normal(), 1.e-7);
    EXPECT_SOFT_NEAR(1.06904497, p.displacement(), 1.e-7);

    // For the normal, the point doesn't matter (we don't check it as that
    // would be expensive)
    auto normal = p.calc_normal(Real3{1.0, 1.0, 1.0});

    EXPECT_SOFT_NEAR(0.26726124, normal[0], 1.e-7);
    EXPECT_SOFT_NEAR(0.53452248, normal[1], 1.e-7);
    EXPECT_SOFT_NEAR(0.80178373, normal[2], 1.e-7);

    cout << p << "\n";
}

//---------------------------------------------------------------------------//

TEST_F(PlaneTest, inline_construction)
{
    real_type d[] = {1.0, 2.0, 3.0, 4.0};

    // Inline construction does not do normalization
    Plane p(d);

    // Check the data
    const real_type* data = p.data();

    EXPECT_EQ(1.0, data[0]);
    EXPECT_EQ(2.0, data[1]);
    EXPECT_EQ(3.0, data[2]);
    EXPECT_EQ(4.0, data[3]);
}

//---------------------------------------------------------------------------//

TEST_F(PlaneTest, translation)
{
    Plane p(Real3(1.0 / sqrt_two, 1.0 / sqrt_two, 0.0), 2 * sqrt_two);
    cout << p << "\n";

    EXPECT_EQ(SignedSense::outside, p.calc_sense(Real3(2.01, 2.01, 0.0)));
    EXPECT_EQ(SignedSense::inside, p.calc_sense(Real3(1.99, 1.99, 0.0)));

    // Translate 2 units to the right
    celeritas::Transform t;
    t.translation(Real3(2, 0, 0));

    auto p2 = p.translated(t);
    cout << p2 << "\n";

    // Check the resulting plane
    EXPECT_VEC_SOFT_EQ(p.normal(), p2.normal());
    EXPECT_EQ(SignedSense::outside, p2.calc_sense(Real3(4.01, 2.01, 0.0)));
    EXPECT_EQ(SignedSense::inside, p2.calc_sense(Real3(3.99, 1.99, 0.0)));
}

//---------------------------------------------------------------------------//

TEST_F(PlaneTest, rotation)
{
    Plane p(Real3(1, 0, 0), 1.25);
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), p.normal());
    EXPECT_SOFT_EQ(1.25, p.displacement());
    cout << p << "\n";

    // 90 degree rotation (quarter turn)
    geometria::Transform t;
    t.rotation(geometria::rotation_matrix(Axis::z, .25));

    // Rotate
    auto p2 = p.transformed(t);
    cout << p2 << "\n";
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), p2.normal());
    EXPECT_SOFT_EQ(1.25, p2.displacement());

    // Rotate again
    auto p3 = p2.transformed(t);
    cout << p3 << "\n";
    EXPECT_VEC_SOFT_EQ(Real3(-1, 0, 0), p3.normal());
    EXPECT_SOFT_EQ(1.25, p3.displacement());
}

//---------------------------------------------------------------------------//
// See notebook gg/python/Surface Tracking Tests.ipynb for references and
// tests
//---------------------------------------------------------------------------//

TEST_F(PlaneTest, Tracking)
{
    // Make a rotated plane in the xy axis
    Plane p(Real3(1 / sqrt_two, 1 / sqrt_two, 0.0), 2 * sqrt_two);

    // Get a point that should have positive sense
    Real3 x({5.41421356, 1.41421356, 0.0});
    EXPECT_EQ(SignedSense::outside, p.calc_sense(x));

    // Calc intersections
    Real3 dir = {-0.70710678, -0.70710678, 0.0};
    normalize_direction(&dir);
    EXPECT_SOFT_NEAR(2.0, calc_intersection(p, x, dir, SurfaceState::off), 1.e-6);

    // Pick a direction such that n\cdot\Omega > 0
    dir = {1.0, 2.0, 3.0};
    normalize_direction(&dir);
    EXPECT_EQ(no_intersection(),
              calc_intersection(p, x, dir, SurfaceState::off));

    // Pick a direction that hits the plane
    dir = {-1, 0.1, 3.0};
    normalize_direction(&dir);
    EXPECT_SOFT_NEAR(9.9430476983098171,
                  calc_intersection(p, x, dir, SurfaceState::off),
                  1.e-6);

    // Place a point on the negative sense
    x = {1.87867966, -2.12132034, 0.0};
    EXPECT_EQ(SignedSense::inside, p.calc_sense(x));

    // Pick a direction such that n\cdot\Omega < 0
    dir = {-1.0, -2.0, 3.0};
    normalize_direction(&dir);
    EXPECT_EQ(no_intersection(),
              calc_intersection(p, x, dir, SurfaceState::off));

    // Pick a direction that hits the plane
    dir = {1, 0.1, 3.0};
    normalize_direction(&dir);
    EXPECT_SOFT_NEAR(12.202831266107504,
                  calc_intersection(p, x, dir, SurfaceState::off),
                  1.e-6);

    // Place a point on the surface
    x = {2.0, 2.0, 0.0};
    EXPECT_EQ(no_intersection(), calc_intersection_from_surface(p, x, dir));
}
