//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstCenteredSphere.cc
 * \brief Tests for class CenteredSphere
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../CenteredSphere.hh"

#include "celeritas_test.hh"
#include "base/ViewField.hh"

using celeritas::CenteredSphere;
using celeritas::no_intersection;
using celeritas::Real3;
using celeritas::SignedSense;
using celeritas::SurfaceState;
using celeritas::Transform;

using cVFDbl = ViewField<const real_type>;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(SphereTest, construction)
{
    EXPECT_EQ(celeritas::SurfaceType::so, CenteredSphere::surface_type());
    EXPECT_EQ(1, CenteredSphere::size());
    EXPECT_EQ(2, CenteredSphere::num_intersections());

    CenteredSphere s(4.0);

    const real_type expected_data[] = {ipow<2>(4)};

    EXPECT_VEC_SOFT_EQ(expected_data, cVFDbl(s.data(), s.data() + s.size()));
}

TEST(SphereTest, translation)
{
    // Set translate-only transform
    Transform t;
    t.translation(Real3(3, -2, 0));

    CenteredSphere s(4.0);

    auto s2 = s.translated(t);
    EXPECT_VEC_SOFT_EQ(Real3(3, -2, 0), s2.origin());
    EXPECT_SOFT_EQ(16.0, s2.radius_sq());
}

TEST(SphereTest, calc_intersections)
{
    real_type dist_iter[CenteredSphere::num_intersections()];

    Real3 trans{-1.1, 2.2, -3.3} CenteredSphere s(4.4);

    // Two intercepts
    {
        Real3 pos{-1.1, -4.4, -3.3};
        Real3 dir{0., 1., 0.};

        pos -= trans;
        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_SOFT_EQ(2.2, dist_iter[0]);
        EXPECT_SOFT_EQ(11.0, dist_iter[1]);
    }
    {
        Real3 pos{10.0, 9.0, 8.0};
        Real3 dir{-1., -1., -1.};

        pos -= trans;
        normalize_direction(&dir);
        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_SOFT_EQ(14.322226209704155, dist_iter[0]);
        EXPECT_SOFT_EQ(19.39502951096999, dist_iter[1]);
    }

    // One intercept
    {
        Real3 pos{-1.1, 2.2, -3.3};
        Real3 dir{0., 0., 1.};

        pos -= trans;
        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_SOFT_EQ(4.4, dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
    {
        Real3 pos{0., 0., 0.};
        Real3 dir{1., 1., 0.};

        pos -= trans;
        normalize_direction(&dir);
        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_SOFT_EQ(2.5170701723978124, dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }

    // No intercepts
    {
        Real3 pos{-5.51, 2.2, -3.3};
        Real3 dir{-1., 0., 0.};

        pos -= trans;
        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_EQ(no_intersection(), dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
    {
        Real3 pos{-1.1, -6.6, 5.51};
        Real3 dir{0., 1., 0.};

        pos -= trans;
        s.calc_intersections(pos, dir, SurfaceState::off, dist_iter);

        EXPECT_EQ(no_intersection(), dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
}

TEST(SphereTest, calc_intersections_on_surface)
{
    real_type dist_iter[CenteredSphere::num_intersections()];

    Real3          trans{4.4, -3.3, 2.2};
    CenteredSphere s(1.1);

    // One intercept
    {
        Real3 pos{3.3, -3.3, 2.2};
        Real3 dir{1., 0., 0.};

        pos -= trans;
        s.calc_intersections(pos, dir, SurfaceState::on, dist_iter);

        EXPECT_SOFT_EQ(2.2, dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
    {
        Real3 pos{4.693987366103667, -2.712025267792666, 3.0819620983110005};
        Real3 dir{-1., -1., -1.};

        pos -= trans;
        normalize_direction(&dir);
        s.calc_intersections(pos, dir, SurfaceState::on, dist_iter);

        EXPECT_SOFT_EQ(2.0368042194996137, dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }

    // No intercept
    {
        Real3 pos{3.3, -3.3, 2.2};
        Real3 dir{-1., 0., 0.};

        pos -= trans;
        s.calc_intersections(pos, dir, SurfaceState::on, dist_iter);

        EXPECT_EQ(no_intersection(), dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
    {
        Real3 pos{4.693987366103667, -2.712025267792666, 3.0819620983110005};
        Real3 dir{1., 1., 1.};

        pos -= trans;
        normalize_direction(&dir);
        s.calc_intersections(pos, dir, SurfaceState::on, dist_iter);

        EXPECT_EQ(no_intersection(), dist_iter[0]);
        EXPECT_EQ(no_intersection(), dist_iter[1]);
    }
}
