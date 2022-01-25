//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstEllipsoidShape.cc
 * \brief Tests for class EllipsoidShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../EllipsoidShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::EllipsoidShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class EllipsoidShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(EllipsoidShapeTest, origin)
{
    EllipsoidShape shape({1, 1, 1});

    EXPECT_EQ("ellipsoid", shape.type());
    EXPECT_SOFT_EQ(4.1887902047863905, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(1, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -1, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1, 1}), this->bbox.upper());

    static const char* expected_surface_names[] = {"sq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 1, 1, 0, 0, 0, -1}),
                       this->get_surface_data(SurfaceType::sq, 0));
}

TEST_F(EllipsoidShapeTest, translated)
{
    EllipsoidShape shape({1, 2, 3});
    EXPECT_SOFT_EQ(3, shape.inradius());

    this->build(shape, Transform{{1, 2, 3}});
    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 4, 6}), this->bbox.upper());

    static const char* expected_surface_names[] = {"sq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({36, 9, 4, -72, -36, -24, 72}),
                       this->get_surface_data(SurfaceType::sq, 0));
}

TEST_F(EllipsoidShapeTest, rotated)
{
    EllipsoidShape shape({3, 2, 1});
    EXPECT_SOFT_EQ(3, shape.inradius());

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-3, -1, -2}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({3, 1, 2}), this->bbox.upper());

    static const char* expected_surface_names[] = {"sq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({4, 36, 9, 0, 0, 0, -36}),
                       this->get_surface_data(SurfaceType::sq, 0));
}
