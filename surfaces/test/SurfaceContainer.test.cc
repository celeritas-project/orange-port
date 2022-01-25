//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstSurfaceContainer.cc
 * \brief SurfaceContainer class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../SurfaceContainer.hh"

#include "celeritas_test.hh"
#include "base/Future.hh"
#include "../Plane.hh"
#include "../PlaneAligned.hh"
#include "../Sphere.hh"
#include "../CylAligned.hh"
#include "../ConeAligned.hh"

using celeritas::SurfaceContainer;
using celeritas::SurfaceType;

//---------------------------------------------------------------------------//
TEST(SurfacesTest, empty)
{
    SurfaceContainer s;
    EXPECT_EQ(0, s.size());
}

//---------------------------------------------------------------------------//
TEST(SurfacesTest, basic)
{
    using namespace celeritas;

    const PlaneX xplane_in(1);
    const PlaneY yplane_in(2);
    const Sphere sphere_in({1, 2, 3}, 4.0);
    const Plane  plane_in({1, 0, 0}, 2.0);

    // Construct surfaces
    SurfaceContainer surfaces;
    surfaces.reserve(4);
    surfaces.push_back(xplane_in);
    surfaces.push_back(yplane_in);
    surfaces.push_back(sphere_in);
    surfaces.push_back(plane_in);

    // Check size/range
    EXPECT_EQ(4, surfaces.size());
    EXPECT_EQ(0, surfaces.all_ids().begin()->get());
    EXPECT_EQ(4, surfaces.all_ids().end()->get());

    // Check types
    EXPECT_EQ(SurfaceType::px, surfaces.get_type(SurfaceId{0}));
    EXPECT_EQ(SurfaceType::py, surfaces.get_type(SurfaceId{1}));
    EXPECT_EQ(SurfaceType::s, surfaces.get_type(SurfaceId{2}));
    EXPECT_EQ(SurfaceType::p, surfaces.get_type(SurfaceId{3}));

    // Check templated types
    ASSERT_TRUE(surfaces.is_type<PlaneX>(SurfaceId{0}));
    ASSERT_TRUE(surfaces.is_type<PlaneY>(SurfaceId{1}));
    ASSERT_TRUE(surfaces.is_type<Sphere>(SurfaceId{2}));
    ASSERT_TRUE(surfaces.is_type<Plane>(SurfaceId{3}));

    // Check that stored data is exactly equal
    {
        auto xplane = surfaces.get<PlaneX>(SurfaceId{0});
        EXPECT_VEC_EQ(xplane_in.view().data, xplane.view().data);

        auto yplane = surfaces.get<PlaneY>(SurfaceId{1});
        EXPECT_VEC_EQ(yplane_in.view().data, yplane.view().data);

        auto sphere = surfaces.get<Sphere>(SurfaceId{2});
        EXPECT_VEC_EQ(sphere_in.view().data, sphere.view().data);

        auto plane = surfaces.get<Plane>(SurfaceId{3});
        EXPECT_VEC_EQ(plane_in.view().data, plane.view().data);
    }

    // Test views
    {
        auto pxview = surfaces.get_view(SurfaceId{0});
        EXPECT_EQ(SurfaceType::px, pxview.type);
        EXPECT_VEC_EQ(xplane_in.view().data, pxview.data);

        auto pview = surfaces.get_view(SurfaceId{3});
        EXPECT_EQ(SurfaceType::p, pview.type);
        EXPECT_VEC_EQ(plane_in.view().data, pview.data);
    }
}

//---------------------------------------------------------------------------//
TEST(SurfacesTest, resize)
{
    using namespace celeritas;

    const PlaneX xplane_in(1);
    const Plane  plane_in({1, 0, 0}, 2.0);

    // Construct surfaces
    SurfaceContainer surfaces;
    surfaces.push_back(xplane_in);
    surfaces.push_back(plane_in);

    // Test resizing
    surfaces.resize(1);
    EXPECT_EQ(1, surfaces.size());
    EXPECT_EQ(SurfaceType::px, surfaces.get_type(SurfaceId{0}));

    // Next placement should be fine
    surfaces.push_back(Sphere({1, 2, 3}, 4.0));
    EXPECT_EQ(SurfaceType::s, surfaces.get_type(SurfaceId{1}));
}
