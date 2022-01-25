//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstExtrudedConvexPolygonShape.cc
 * \brief Tests for class ExtrudedConvexPolygonShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../ExtrudedConvexPolygonShape.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::ExtrudedConvexPolygonShape;

using constants::pi;
using constants::sqrt_two;
using geometria::rotation_matrix;
using std::cos;
using std::tan;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class ExtrudedConvexPolygonShapeTest : public ::celeritas::ShapeTest
{
  protected:
    //// TYPEDEFS ////
    using Point    = Array<real_type, 2>;
    using VecPoint = std::vector<Point>;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(ExtrudedConvexPolygonShapeTest, triangle)
{
    const real_type apothem  = 1.0;
    const real_type r        = apothem / cos(pi / 3);
    const real_type x        = r * cos(pi / 6);
    VecPoint        vertices = {{-x, -apothem}, {x, -apothem}, {0, r}};
    ExtrudedConvexPolygonShape shape(vertices, -1.2, 2.4);

    this->build(shape);

    static const logic_int expected_logic[] = {0l,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1.73205080756888, -1, -1.2}),
                       this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1.73205080756888, 2, 2.4}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::py, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, 0.5, 0, 1}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, -1}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1.2}),
                       this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({2.4}),
                       this->get_surface_data(SurfaceType::pz, 4));
}

TEST_F(ExtrudedConvexPolygonShapeTest, translated)
{
    const real_type apothem  = 1.0;
    const real_type r        = apothem / cos(pi / 3);
    const real_type x        = r * cos(pi / 6);
    VecPoint        vertices = {{-x, -apothem}, {x, -apothem}, {0, r}};
    ExtrudedConvexPolygonShape shape(vertices, -1.2, 2.4);

    this->build(shape, Transform{{1, 2, 3}});

    static const logic_int expected_logic[] = {0l,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-0.732050807568877, 1, 1.8}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2.73205080756888, 4, 5.4}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::py, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, 0.5, 0, 2.8660254037844}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, -1.1339745962156}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({1.8}),
                       this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({5.4}),
                       this->get_surface_data(SurfaceType::pz, 4));
}

TEST_F(ExtrudedConvexPolygonShapeTest, rotated)
{
    const real_type apothem  = 1.0;
    const real_type r        = apothem / cos(pi / 3);
    const real_type x        = r * cos(pi / 6);
    VecPoint        vertices = {{-x, -apothem}, {x, -apothem}, {0, r}};
    ExtrudedConvexPolygonShape shape(vertices, -1.2, 2.4);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});

    static const logic_int expected_logic[] = {0l,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1.7320508075689, -2.4, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1.7320508075689, 1.2, 2}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, 0, 0.5, 1}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, 0, -0.5, -1}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({1.2}),
                       this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2.4}),
                       this->get_surface_data(SurfaceType::py, 4));
}

TEST_F(ExtrudedConvexPolygonShapeTest, square)
{
    real_type                  half_width = 1.0;
    VecPoint                   vertices   = {{-half_width, -half_width},
                         {half_width, -half_width},
                         {half_width, half_width},
                         {-half_width, half_width}};
    ExtrudedConvexPolygonShape shape(vertices, -1.23, 2.34);

    this->build(shape);

    EXPECT_VEC_EQ(CSGCell::from_string("0 1 ~ & 2 ~ & 3 & 4 & 5 ~ & ").logic(),
                  this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -1, -1.23}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1, 2.34}), this->bbox.upper());

    static const char* const expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::py, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::px, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1.23}),
                       this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({2.34}),
                       this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ExtrudedConvexPolygonShapeTest, hexagon)
{
    constexpr static real_type half_sqrt_three = 0.8660254037844386;

    const real_type apothem = 2;
    const real_type r       = apothem / cos(pi / 6);
    const real_type y       = r * sin(pi / 6);

    VecPoint                   vertices = {{apothem, -y},
                         {apothem, y},
                         {0, r},
                         {-apothem, y},
                         {-apothem, -y},
                         {0, -r}};
    ExtrudedConvexPolygonShape shape(vertices, -1, 1);

    this->build(shape);

    EXPECT_VEC_EQ(
        CSGCell::from_string("0 ~ 1 ~ & 2 & 3 & 4 & 5 ~ & 6 & 7 ~ & ").logic(),
        this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -2.3094010767585, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 2.3094010767585, 1}), this->bbox.upper());

    static const char* const expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "p4", "p5", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(8, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.86602540378444, 0, 2}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.86602540378444, 0, -2}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2}),
                       this->get_surface_data(SurfaceType::px, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.86602540378444, 0, -2}),
                       this->get_surface_data(SurfaceType::p, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.86602540378444, 0, 2}),
                       this->get_surface_data(SurfaceType::p, 5));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 6));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 7));
}

TEST_F(ExtrudedConvexPolygonShapeTest, generic_polygon)
{
    VecPoint vertices = {{-1, 1}, {-1, -1}, {1, 0}, {0.5, 1.5}};
    ExtrudedConvexPolygonShape shape(vertices, -1, 1);

    this->build(shape);

    EXPECT_VEC_EQ(CSGCell::from_string("0 1 ~ & 2 ~ & 3 & 4 & 5 ~ & ").logic(),
                  this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -1, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1.5, 1}), this->bbox.upper());

    static const char* const expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.44721359549996, -0.89442719099992, 0, 0.44721359549996}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.94868329805051, 0.31622776601684, 0, 0.94868329805051}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.31622776601684, -0.94868329805051, 0, -1.2649110640674}),
        this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 5));
}
