//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstCylinderSegmentShape.cc
 * \brief Tests for class CylinderSegmentShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../CylinderSegmentShape.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::CylinderSegmentShape;

using constants::pi;
using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class CylinderSegmentShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CylinderSegmentShapeTest, origin)
{
    CylinderSegmentShape shape(1.1, 2.2, 0, 1.5708, -1.0, 1.0);

    EXPECT_EQ("cylinder_segment", shape.type());
    EXPECT_SOFT_EQ(5.702004, shape.volume());
    EXPECT_FALSE(shape.is_convex());
    EXPECT_SOFT_EQ(3.2680269276736387, shape.inradius());

    this->build(shape);

    // -3 = LOGICAL_AND, -1 = LOGICAL_NOT
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
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2.2, 2.2, 1}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"ciz", "coz", "mz", "pz", "p0", "p1"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1.21}),
                       this->get_surface_data(SurfaceType::czo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({4.84}),
                       this->get_surface_data(SurfaceType::czo, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::px, 5));
}

TEST_F(CylinderSegmentShapeTest, translated)
{
    CylinderSegmentShape shape(0.1, 2.2, pi / 2, pi, -1.1, 1.3);

    EXPECT_SOFT_EQ(3.3346664001066135, shape.inradius());

    this->build(shape, Transform{{1, 2, 3}});

    // -3 = LOGICAL_AND, -1 = LOGICAL_NOT
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
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1.2, -0.20000000000002, 1.9}),
                       this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 4.2, 4.3}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"ciz", "coz", "mz", "pz", "p0", "p1"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 2, 0.01}),
                       this->get_surface_data(SurfaceType::cz, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 2, 4.84}),
                       this->get_surface_data(SurfaceType::cz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1.9}),
                       this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({4.3}),
                       this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::px, 4));
}

TEST_F(CylinderSegmentShapeTest, rotated)
{
    CylinderSegmentShape shape(0.1, 2.2, pi, pi / 4, 0, 1.3);
    EXPECT_SOFT_EQ(3.1784430150625638, shape.inradius());

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});

    // -3 = LOGICAL_AND, -1 = LOGICAL_NOT
    static const logic_int expected_logic[] = {0l,
                                               1l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2.2, -1.3, -2.2}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2.2, 0, 0}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"ciz", "coz", "mz", "pz", "p0", "p1"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.01}),
                       this->get_surface_data(SurfaceType::cyo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({4.84}),
                       this->get_surface_data(SurfaceType::cyo, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1.3}),
                       this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.70710678118655, 0, -0.70710678118655, 0}),
                       this->get_surface_data(SurfaceType::p, 5));
}
