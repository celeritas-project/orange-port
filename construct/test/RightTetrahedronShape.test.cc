//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstRightTetrahedronShape.cc
 * \brief Tests for class RightTetrahedronShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../RightTetrahedronShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::RightTetrahedronShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class RightTetrahedronShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(RightTetrahedronShapeTest, origin)
{
    RightTetrahedronShape shape(1, 1, 1);

    EXPECT_EQ("right_tetrahedron", shape.type());
    EXPECT_SOFT_EQ(0.16666666666666666, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[]
        = {0l, 1l, LOGIC_AND, 2l, LOGIC_AND, 3l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1, 1}), this->bbox.upper());

    static const char* expected_surface_names[] = {"mx", "my", "mz", "p0"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(4, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.57735026918963,
                               0.57735026918963,
                               0.57735026918963,
                               0.57735026918963}),
                       this->get_surface_data(SurfaceType::p, 3));
}

TEST_F(RightTetrahedronShapeTest, offdiagonal)
{
    RightTetrahedronShape shape(1, 2, 3);

    EXPECT_EQ("right_tetrahedron", shape.type());
    EXPECT_SOFT_EQ(1, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[]
        = {0l, 1l, LOGIC_AND, 2l, LOGIC_AND, 3l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 2, 3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"mx", "my", "mz", "p0"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(4, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.85714285714286,
                               0.42857142857143,
                               0.28571428571429,
                               0.85714285714286}),
                       this->get_surface_data(SurfaceType::p, 3));
}

TEST_F(RightTetrahedronShapeTest, translated)
{
    RightTetrahedronShape shape(1, 2, 3);

    this->build(shape, Transform{{1, 2, 3}});
    static const logic_int expected_logic[]
        = {0l, 1l, LOGIC_AND, 2l, LOGIC_AND, 3l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({1, 2, 3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 4, 6}), this->bbox.upper());

    static const char* expected_surface_names[] = {"mx", "my", "mz", "p0"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(4, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.85714285714286,
                               0.42857142857143,
                               0.28571428571429,
                               3.4285714285714}),
                       this->get_surface_data(SurfaceType::p, 3));
}

TEST_F(RightTetrahedronShapeTest, rotated)
{
    RightTetrahedronShape shape(1, 2, 3);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
    static const logic_int expected_logic[] = {
        0l, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND, 3l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, -3, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 0, 2}), this->bbox.upper());

    static const char* expected_surface_names[] = {"mx", "my", "mz", "p0"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(4, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.85714285714286,
                               -0.28571428571429,
                               0.42857142857143,
                               0.85714285714286}),
                       this->get_surface_data(SurfaceType::p, 3));
}
