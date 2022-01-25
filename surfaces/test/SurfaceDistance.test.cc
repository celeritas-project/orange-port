//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstSurfaceDistance.cc
 * \brief Tests for class SurfaceDistance
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../SurfaceDistance.hh"

#include "celeritas_test.hh"
#include "../PlaneAligned.hh"
#include "../Plane.hh"
#include "../CylCentered.hh"
#include "../CenteredSphere.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(SurfaceDistanceTest, ortho_plane)
{
    using Dist_t = SurfaceDistance<PlaneX>;

    PlaneX p1(10.0);
    PlaneX p2(10.0 + 1e-3);

    EXPECT_SOFT_EQ(1e-3, Dist_t::global(p1, p2));
    EXPECT_SOFT_EQ(1e-3, Dist_t::global(p2, p1));
}

//---------------------------------------------------------------------------//

TEST(SurfaceDistanceTest, plane)
{
    using Dist_t = SurfaceDistance<Plane>;

    Real3     n(1.0, 2.0, 3.0);
    real_type d = 4.0 / vector_magnitude(n);
    normalize_direction(&n);

    Plane p1(n, d);
    Plane p2(n, d + 1e-3);

    EXPECT_SOFT_EQ(1e-3, Dist_t::global(p1, p2));
    EXPECT_SOFT_EQ(1e-3, Dist_t::global(p2, p1));
}

//---------------------------------------------------------------------------//

TEST(SurfaceDistanceTest, cyl)
{
    using Cyl_t  = CylCentered<Axis::x>;
    using Dist_t = SurfaceDistance<Cyl_t>;

    Cyl_t c1(2.0);
    Cyl_t c2(2.0 + 1e-3);

    EXPECT_SOFT_EQ(1e-3, Dist_t::global(c1, c2));
    EXPECT_SOFT_EQ(1e-3, Dist_t::global(c2, c1));
}
