//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstPrismShape.cc
 * \brief Tests for class PrismShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../PrismShape.hh"

#include <cmath>

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::PrismShape;

using geometria::rotation_matrix;

using constants::pi;
using constants::sqrt_two;
using std::cos;

constexpr static real_type half_sqrt_three = 0.8660254037844386;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class PrismShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(PrismShapeTest, triangle)
{
    real_type  apothem = 1.0;
    PrismShape shape(3, apothem, 0.0, -1.2, 2.4);

    EXPECT_EQ("prism", shape.type());
    EXPECT_SOFT_EQ(18.706148721743865, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(1.0, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {0l,
                                               LOGIC_NOT,
                                               1l,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -1, -1.2}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 2, 2.4}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, 0.5, 0, 1}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, -1}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1.2}),
                       this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({2.4}),
                       this->get_surface_data(SurfaceType::pz, 4));
}

TEST_F(PrismShapeTest, square)
{
    real_type  half_width = 1.0;
    PrismShape shape(4, half_width, 0.0, -1.23, 2.34);

    EXPECT_EQ("prism", shape.type());
    EXPECT_SOFT_EQ(14.279999999999998, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(1, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {0l,
                                               LOGIC_NOT,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -1, -1.23}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1, 2.34}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::px, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1.23}),
                       this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({2.34}),
                       this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(PrismShapeTest, pentagon)
{
    real_type  apothem = 1.5;
    PrismShape shape(5, apothem, 0.0, -1, 1);

    EXPECT_EQ("prism", shape.type());
    EXPECT_SOFT_EQ(16.34720688012062, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(1.5, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {0l,
                                               LOGIC_NOT,
                                               1l,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_AND,
                                               6l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1.8541019662497, -1.5, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1.8541019662497, 1.8541019662497, 1}),
                       this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "p4", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(7, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.58778525229247, 0.80901699437495, 0, 1.5}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.58778525229247, -0.80901699437495, 0, -1.5}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.95105651629515, 0.30901699437495, 0, -1.5}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1.5}),
                       this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.95105651629515, -0.30901699437495, 0, 1.5}),
                       this->get_surface_data(SurfaceType::p, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 5));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 6));
}

// Test the KENO equivalent of a RHEXPRISM
TEST_F(PrismShapeTest, hexagon)
{
    real_type  apothem = 2;
    PrismShape shape(6, apothem, 0.0, -1, 1);

    this->build(shape);
    static const logic_int expected_logic[] = {0l,
                                               LOGIC_NOT,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               6l,
                                               LOGIC_AND,
                                               7l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2.3094010767585, -2, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2.3094010767585, 2, 1}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "p4", "p5", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(8, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, 0.5, 0, 2}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, -2}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, 0.5, 0, -2}),
                       this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2}),
                       this->get_surface_data(SurfaceType::py, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, 2}),
                       this->get_surface_data(SurfaceType::p, 5));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 6));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 7));
}

// Test the KENO equivalent of a HEXPRISM (pointy-top)
TEST_F(PrismShapeTest, rot_hexagon)
{
    real_type  apothem = 2;
    PrismShape shape(6, apothem, 0.5, -1, 1);

    EXPECT_EQ("prism", shape.type());
    EXPECT_SOFT_EQ(27.712812921102028, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(2, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {0l,
                                               LOGIC_NOT,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               6l,
                                               LOGIC_AND,
                                               7l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -2.3094010767585, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 2.3094010767585, 1}), this->bbox.upper());

    static const char* expected_surface_names[]
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

TEST_F(PrismShapeTest, translated)
{
    real_type  apothem = 2;
    PrismShape shape(6, apothem, 0.5, -1, 1);

    this->build(shape, Transform{{1, 2, 3}});
    static const logic_int expected_logic[] = {0l,
                                               LOGIC_NOT,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               6l,
                                               LOGIC_AND,
                                               7l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -0.3094010767585, 2}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({3, 4.3094010767585, 4}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "p4", "p5", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(8, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.86602540378444, 0, 4.2320508075689}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.86602540378444, 0, -3.2320508075689}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::px, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.86602540378444, 0, 0.23205080756888}),
                       this->get_surface_data(SurfaceType::p, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.86602540378444, 0, 0.76794919243112}),
                       this->get_surface_data(SurfaceType::p, 5));
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::pz, 6));
    EXPECT_VEC_SOFT_EQ(VecDbl({4}), this->get_surface_data(SurfaceType::pz, 7));
}

TEST_F(PrismShapeTest, rotated)
{
    real_type  apothem = 2;
    PrismShape shape(6, apothem, 0.5, -1, 1);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
    static const logic_int expected_logic[] = {0l,
                                               LOGIC_NOT,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               6l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               7l,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -1, -2.3094010767585}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 1, 2.3094010767585}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "p4", "p5", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(8, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0, 0.86602540378444, 2}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0, -0.86602540378444, -2}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2}),
                       this->get_surface_data(SurfaceType::px, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0, 0.86602540378444, -2}),
                       this->get_surface_data(SurfaceType::p, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0, -0.86602540378444, 2}),
                       this->get_surface_data(SurfaceType::p, 5));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::py, 6));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::py, 7));
}
