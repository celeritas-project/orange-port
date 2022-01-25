//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstCylinderShape.cc
 * \brief Tests for class CylinderShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../CylinderShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "base/Definitions.hh"
#include "orange/surfaces/CylAligned.hh"
#include "orange/surfaces/CylCentered.hh"
#include "ShapeTest.hh"

using celeritas::CylinderShape;

using celeritas::CylAligned;
using celeritas::CylCentered;
using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class CylinderShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CylinderShapeTest, origin)
{
    CylinderShape shape(Axis::z, 1, -1, 1);

    EXPECT_EQ("cylinder", shape.type());
    EXPECT_SOFT_EQ(6.2831853071795862, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(1.0, shape.inradius());

    EXPECT_EQ(Axis::z, shape.axis());
    EXPECT_EQ(-1, shape.lo());
    EXPECT_EQ(1, shape.hi());
    EXPECT_EQ(1, shape.radius());

    this->build(shape);

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -1, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1, 1}), this->bbox.upper());

    static const char* expected_surface_names[] = {"coz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}),
                       this->get_surface_data(SurfaceType::czo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 2));
}

TEST_F(CylinderShapeTest, origin_z_aligned)
{
    CylinderShape shape(1, 2);

    EXPECT_EQ("cylinder", shape.type());
    EXPECT_SOFT_EQ(6.2831853071795862, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(1.0, shape.inradius());

    EXPECT_EQ(Axis::z, shape.axis());
    EXPECT_EQ(-1, shape.lo());
    EXPECT_EQ(1, shape.hi());
    EXPECT_EQ(1, shape.radius());

    this->build(shape);

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -1, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1, 1}), this->bbox.upper());

    static const char* expected_surface_names[] = {"coz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}),
                       this->get_surface_data(SurfaceType::czo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 2));
}

TEST_F(CylinderShapeTest, origin_z_aligned_inf)
{
    CylinderShape shape(2, inf);

    EXPECT_EQ("cylinder", shape.type());
    EXPECT_SOFT_EQ(inf, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(2.0, shape.inradius());

    EXPECT_EQ(Axis::z, shape.axis());
    EXPECT_SOFT_EQ(-inf, shape.lo());
    EXPECT_SOFT_EQ(inf, shape.hi());
    EXPECT_SOFT_EQ(2.0, shape.radius());

    this->build(shape);

    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -2, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 2, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"coz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({4}),
                       this->get_surface_data(SurfaceType::czo, 0));
}

TEST_F(CylinderShapeTest, translated)
{
    CylinderShape shape(Axis::x, 1.3, -1.2, 2.1);

    this->build(shape, Transform{{1, 2, 3}});

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-0.2, 0.69999999999999, 1.7}),
                       this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({3.1, 3.3, 4.3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"cox", "mx", "px"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({2, 3, 1.69}),
                       this->get_surface_data(SurfaceType::cx, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-0.2}),
                       this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({3.1}),
                       this->get_surface_data(SurfaceType::px, 2));
}

TEST_F(CylinderShapeTest, rotated)
{
    CylinderShape shape(Axis::y, 3.1, -2.1, 1.2);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-3.1, -3.1, -2.1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({3.1, 3.1, 1.2}), this->bbox.upper());

    static const char* expected_surface_names[] = {"coy", "my", "py"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({9.61}),
                       this->get_surface_data(SurfaceType::czo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2.1}),
                       this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1.2}),
                       this->get_surface_data(SurfaceType::pz, 2));
}
