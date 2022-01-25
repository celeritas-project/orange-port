//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstSurfaceAction.cc
 * \brief Tests for class SurfaceAction
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../SurfaceAction.hh"

#include "celeritas_test.hh"

using celeritas::make_surface_action;
using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class SurfaceActionTest : public ::Test
{
  protected:
    //// TYPE ALIASES ////

  protected:
    void SetUp()
    {
        surfaces.push_back(PlaneX(1));
        surfaces.push_back(PlaneY(2));
        surfaces.push_back(PlaneZ(3));
        surfaces.push_back(CenteredSphere(4));
        surfaces.push_back(CCylX(5));
        surfaces.push_back(CCylY(6));
        surfaces.push_back(CCylZ(7));
        surfaces.push_back(Plane({1, 1, 0}, {1, 0, 0}));
        surfaces.push_back(Sphere({1, 2, 3}, 4));
        surfaces.push_back(CylX({1, 2, 3}, 5));
        surfaces.push_back(CylY({1, 2, 3}, 6));
        surfaces.push_back(CylZ({1, 2, 3}, 7));
        surfaces.push_back(ConeX({2, 3, 4}, .05));
        surfaces.push_back(ConeY({2, 3, 4}, .04));
        surfaces.push_back(ConeZ({2, 3, 4}, .03));
        surfaces.push_back(
            SimpleQuadric({1, 2, 3}, {4, 5, 6}, 7, {-1, -2, -3}));
        surfaces.push_back(GeneralQuadric({0, 1, 2}, {3, 4, 5}, {6, 7, 8}, 9));
    }

  protected:
    //// DATA ////
    SurfaceContainer surfaces;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

struct ToString
{
    template<class S>
    std::string operator()(S&& surf) const
    {
        std::ostringstream os;
        os << surf;
        return os.str();
    }
};

TEST_F(SurfaceActionTest, test_strings)
{
    // Create functor
    auto surf_to_string = make_surface_action(surfaces, ToString());

    // Loop over all surfaces and apply
    std::vector<std::string> strings;
    for (auto id : surfaces.all_ids())
    {
        strings.push_back(surf_to_string(id));
    }

    // clang-format off
    const std::string expected_strings[] = {
        "Plane: x=1",
        "Plane: y=2",
        "Plane: z=3",
        "Sphere: r=4",
        "Cyl x: r=5",
        "Cyl y: r=6",
        "Cyl z: r=7",
        "Plane: n=(0.707107 0.707107 0), "
        "d=0.707107",
        "Sphere: r=4 at 1 2 3",
        "Cyl x: r=5, y=2, z=3",
        "Cyl y: r=6, x=1, z=3",
        "Cyl z: r=7, x=1, y=2",
        "Cone x: tangent=0.05 at 2 3 4",
        "Cone y: tangent=0.04 at 2 3 4",
        "Cone z: tangent=0.03 at 2 3 4",
        "SQuadric: {1,2,3} {6,13,24} 107",
        "GQuadric: {0,1,2} {3,4,5} {6,7,8} 9"};
    // clang-format on
    EXPECT_VEC_EQ(expected_strings, strings);
}

//---------------------------------------------------------------------------//

struct AccumulateIntersections
{
    int* num_intersections;

    template<class S>
    void operator()(const S& surf)
    {
        *num_intersections += S::num_intersections();
    }
};

TEST_F(SurfaceActionTest, test_intersections)
{
    // Create functor
    int  num_intersections  = 0;
    auto accum_intersection = make_surface_action(
        surfaces, AccumulateIntersections{&num_intersections});

    // Loop over all surfaces and apply
    for (auto id : surfaces.all_ids())
    {
        accum_intersection(id);
    }
    EXPECT_EQ(30, num_intersections);
}
