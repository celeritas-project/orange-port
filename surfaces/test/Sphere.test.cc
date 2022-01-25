//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstSphere.cc
 * \brief Sphere class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Sphere.hh"

#include "celeritas_test.hh"
#include "base/ViewField.hh"
#include "orange/Transform.hh"

using celeritas::Sphere;

using celeritas::no_intersection;
using celeritas::Real3;
using celeritas::SignedSense;
using celeritas::SurfaceState;
using celeritas::Transform;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(SphereTest, construction)
{
    EXPECT_EQ(celeritas::SurfaceType::s, Sphere::surface_type());
    EXPECT_EQ(4, Sphere::size());
    EXPECT_EQ(2, Sphere::num_intersections());

    Sphere s(Real3(1, 2, 3), 4.0);

    const real_type expected_data[] = {1, 2, 3, ipow<2>(4)};

    EXPECT_VEC_SOFT_EQ(
        expected_data,
        ViewConstField<real_type>(s.data(), s.data() + s.size()));
}

TEST(SphereTest, translation)
{
    // Set translate-only transform
    Transform t;
    t.translation(Real3(3, -2, 0));

    Sphere s(Real3(1, 2, 3), 4.0);

    auto s2 = s.translated(t);
    EXPECT_VEC_SOFT_EQ(Real3(4, 0, 3), s2.origin());
    EXPECT_SOFT_EQ(16.0, s2.radius_sq());
}

TEST(SphereTest, calc_intersections)
{
    real_type dist_iter[Sphere::num_intersections()];

    Sphere s(Real3{-1.1, 2.2, -3.3}, 4.4);

    // Two intercepts
    {
        Real3 pos{-1.1, -4.4, -3.3};
        Real3 dir{0., 1., 0.};

        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_SOFT_EQ(2.2, dist_iter[0]);
        EXPECT_SOFT_EQ(11.0, dist_iter[1]);
    }
    {
        Real3 pos{10.0, 9.0, 8.0};
        Real3 dir{-1., -1., -1.};

        normalize_direction(&dir);
        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_SOFT_EQ(14.322226209704155, dist_iter[0]);
        EXPECT_SOFT_EQ(19.39502951096999, dist_iter[1]);
    }

    // One intercept
    {
        Real3 pos{-1.1, 2.2, -3.3};
        Real3 dir{0., 0., 1.};

        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_SOFT_EQ(4.4, dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
    {
        Real3 pos{0., 0., 0.};
        Real3 dir{1., 1., 0.};

        normalize_direction(&dir);
        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_SOFT_EQ(2.5170701723978124, dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }

    // No intercepts
    {
        Real3 pos{-5.51, 2.2, -3.3};
        Real3 dir{-1., 0., 0.};

        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_EQ(no_intersection(), dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
    {
        Real3 pos{-1.1, -6.6, 5.51};
        Real3 dir{0., 1., 0.};

        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_EQ(no_intersection(), dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
}

TEST(SphereTest, calc_intersections_on_surface)
{
    real_type dist_iter[Sphere::num_intersections()];

    Sphere s(Real3{4.4, -3.3, 2.2}, 1.1);

    // One intercept
    {
        Real3 pos{3.3, -3.3, 2.2};
        Real3 dir{1., 0., 0.};

        s.calc_intersections(pos, dir, SurfaceState::on, dist_iter);

        EXPECT_SOFT_EQ(2.2, dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
    {
        Real3 pos{4.693987366103667, -2.712025267792666, 3.0819620983110005};
        Real3 dir{-1., -1., -1.};

        normalize_direction(&dir);
        s.calc_intersections(pos, dir, SurfaceState::on, dist_iter);

        EXPECT_SOFT_EQ(2.0368042194996137, dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }

    // No intercept
    {
        Real3 pos{3.3, -3.3, 2.2};
        Real3 dir{-1., 0., 0.};

        s.calc_intersections(pos, dir, SurfaceState::on, dist_iter);

        EXPECT_EQ(no_intersection(), dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
    {
        Real3 pos{4.693987366103667, -2.712025267792666, 3.0819620983110005};
        Real3 dir{1., 1., 1.};

        normalize_direction(&dir);
        s.calc_intersections(pos, dir, SurfaceState::on, dist_iter);

        EXPECT_EQ(no_intersection(), dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
}

TEST(SphereTest, calc_intersections_degenerate)
{
    real_type distances[Sphere::num_intersections()];

    const real_type eps = 1.e-4;
    Sphere          s(Real3{1.0, 0.0, 0}, 1.0);

    // Heading toward
    s.calc_intersections(
        Real3{eps, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);
    EXPECT_SOFT_NEAR(2.0 - eps, distances[0], eps);
    EXPECT_EQ(no_intersection(), distances[1]);

    // Heading away, slightly inside
    s.calc_intersections(
        Real3{eps, 0, 0}, Real3{-1, 0, 0}, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(eps, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // Tangent
    s.calc_intersections(
        Real3{eps, 0, 0}, Real3{0, 1, 0}, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // Heading in
    s.calc_intersections(
        Real3{eps, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);
    EXPECT_SOFT_NEAR(2.0 + eps, distances[0], eps);
    EXPECT_EQ(no_intersection(), distances[1]);

    // Heading away
    s.calc_intersections(
        Real3{2.0 - eps, 0, 0}, Real3{1, 0, 0}, SurfaceState::off, distances);
    EXPECT_SOFT_NEAR(eps, distances[0], eps);
    EXPECT_EQ(no_intersection(), distances[1]);

    // Tangent
    s.calc_intersections(
        Real3{eps, 0, 0}, Real3{0, 1, 0}, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(std::sqrt(2 * eps * 1.0 - ipow<2>(eps)), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

//---------------------------------------------------------------------------//
/*!
 * This is a failure point of the "coincident" point in keno_spheres
 */
TEST(TrickySphereTest, top_downward)
{
    const Sphere left(Real3(-4.5, -4.5, 0), 10.0);
    const Sphere right(Real3(4.5, -4.5, 0), 10.0);
    const Real3  dir(0, -1, 0);

    Real3     pos;
    real_type distances[Sphere::num_intersections()];

    // Initial position
    pos = {0, 10, 0};
    EXPECT_EQ(SignedSense::outside, left.calc_sense(pos));
    EXPECT_EQ(SignedSense::outside, right.calc_sense(pos));

    // Initial distances
    real_type expected_distance = 5.5697144502541249;
    left.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(expected_distance, distances[0]);
    EXPECT_SOFT_EQ(23.430285549745875, distances[1]);

    right.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(expected_distance, distances[0]);
    EXPECT_SOFT_EQ(23.430285549745875, distances[1]);

    // It crosses outside the left sphere
    pos += dir * expected_distance;
    expected_distance = 17.86057109949175;
    left.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(expected_distance, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // We'd expect to be "on" the second sphere but floating point error puts
    // us just inside it.
    EXPECT_EQ(SignedSense::inside, right.calc_sense(pos));
    right.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(0.0, distances[0]);
    EXPECT_SOFT_EQ(expected_distance, distances[1]);

    // Transporting the distance of zero puts us "on" the surface
    right.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(expected_distance, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // Transport to the other side of the sphere
    pos += dir * distances[0];
    right.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

/*!
 * This is part A of the "coincident" point in keno_spheres giving
 * problems.
 */
TEST(TrickySphereTest, bottom_downward)
{
    const Sphere left(Real3(-4.5, -4.5, 0), 10.0);
    const Sphere right(Real3(4.5, -4.5, 0), 10.0);
    const Real3  dir(0, -1, 0);

    Real3     pos;
    real_type distances[Sphere::num_intersections()];

    // Initial position
    pos = {0, -5, 0};
    EXPECT_EQ(SignedSense::inside, left.calc_sense(pos));
    EXPECT_EQ(SignedSense::inside, right.calc_sense(pos));

    // Initial distances
    real_type expected_distance = 8.4302855497458751;
    left.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(expected_distance, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    right.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(expected_distance, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // It crosses outside the left sphere
    pos[Axis::y] += -1 * expected_distance;
    left.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_EQ(no_intersection(), distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    // We'd expect to be "on" the second sphere but floating point error puts
    // us just inside it
    // EXPECT_EQ(SignedSense::outside, right.calc_sense(pos));
    EXPECT_EQ(SignedSense::inside, right.calc_sense(pos));
    right.calc_intersections(pos, dir, SurfaceState::off, distances);
    // Because we return "no intersection" instead of a small positive
    // distance, a particle gets lost at this sphere boundary
    // EXPECT_SOFT_NEAR(1.e-14, distances[0], 0.5);
    EXPECT_EQ(no_intersection(), distances[1]);
    EXPECT_EQ(no_intersection(), distances[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Test two exactly tangent spheres (at the origin) inside a global
 * sphere
 * 0: -A
 * 1: -B
 * 2: +A +B -C
 * EXT: +C
 */
TEST(TangentSpheres, left_rightward)
{
    const Sphere left({-1, 0, 0}, 1.0);
    const Sphere right({1, 0, 0}, 1.0);
    const Sphere global({0, 0, 0}, 3.0);

    Real3     pos;
    Real3     dir{1, 0, 0};
    real_type distances[Sphere::num_intersections()];

    // Initial position
    pos = {-3, 0, 0};
    EXPECT_EQ(SignedSense::outside, left.calc_sense(pos));
    EXPECT_EQ(SignedSense::outside, right.calc_sense(pos));
    EXPECT_EQ(SignedSense::on, global.calc_sense(pos));

    global.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(0.0, distances[0]);
    EXPECT_EQ(6.0, distances[1]);

    // Move: distance of zero, now 'on' global surface
    global.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(6.0, distances[0]);
    EXPECT_EQ(no_intersection(), distances[1]);

    left.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(1.0, distances[0]);
    EXPECT_SOFT_EQ(3.0, distances[1]);

    // Move to left side of left sphere
    pos[Axis::x] += distances[0];
    EXPECT_EQ(SignedSense::on, left.calc_sense(pos));
    EXPECT_EQ(SignedSense::outside, right.calc_sense(pos));
    EXPECT_EQ(SignedSense::inside, global.calc_sense(pos));

    left.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(2.0, distances[0]);
    EXPECT_SOFT_EQ(no_intersection(), distances[1]);

    // Move to right side of left sphere
    pos[Axis::x] += distances[0];
    EXPECT_SOFT_EQ(0.0, pos[Axis::x]);
    EXPECT_EQ(SignedSense::on, left.calc_sense(pos));
    EXPECT_EQ(SignedSense::on, right.calc_sense(pos));
    EXPECT_EQ(SignedSense::inside, global.calc_sense(pos));
    left.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(no_intersection(), distances[0]);
    EXPECT_SOFT_EQ(no_intersection(), distances[1]);

    // Crossed surface: in adjacent 'global' cell (pos pos neg -> inside)
    left.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(no_intersection(), distances[0]);
    right.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(0.0, distances[0]);
    global.calc_intersections(pos, dir, SurfaceState::off, distances);
    EXPECT_SOFT_EQ(3.0, distances[0]);

    // Move distance of zero to right surface, cross inside.
    right.calc_intersections(pos, dir, SurfaceState::on, distances);
    EXPECT_SOFT_EQ(2.0, distances[0]);
}
