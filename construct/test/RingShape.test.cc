//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstRingShape.cc
 * \brief Tests for class RingShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../RingShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::RingShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class RingShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(RingShapeTest, origin)
{
    RingShape shape(1.0, 4.0, -2.0, 2.0);

    EXPECT_EQ("ring", shape.type());
    EXPECT_SOFT_EQ(188.49555921538757, shape.volume());
    EXPECT_FALSE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {
        0l, 1l, LOGIC_NOT, LOGIC_AND, 2l, LOGIC_AND, 3l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-4, -4, -2}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4, 4, 2}), this->bbox.upper());

    static const char* expected_surface_names[] = {"ciz", "coz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(4, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}),
                       this->get_surface_data(SurfaceType::czo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({16}),
                       this->get_surface_data(SurfaceType::czo, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2}),
                       this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::pz, 3));
}

TEST_F(RingShapeTest, origin_simple)
{
    RingShape shape(1.0, 4.0, 4.0);

    EXPECT_EQ("ring", shape.type());
    EXPECT_SOFT_EQ(188.49555921538757, shape.volume());
    EXPECT_FALSE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {
        0l, 1l, LOGIC_NOT, LOGIC_AND, 2l, LOGIC_AND, 3l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-4, -4, -2}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4, 4, 2}), this->bbox.upper());

    static const char* expected_surface_names[] = {"ciz", "coz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(4, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}),
                       this->get_surface_data(SurfaceType::czo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({16}),
                       this->get_surface_data(SurfaceType::czo, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2}),
                       this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::pz, 3));
}

TEST_F(RingShapeTest, translated)
{
    RingShape shape(1.0, 4.0, 0.0, 2.0);

    this->build(shape, Transform{{1, 2, 3}});
    static const logic_int expected_logic[] = {
        0l, 1l, LOGIC_NOT, LOGIC_AND, 2l, LOGIC_AND, 3l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-3, -2, 3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({5, 6, 5}), this->bbox.upper());

    static const char* expected_surface_names[] = {"ciz", "coz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(4, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 2, 1}),
                       this->get_surface_data(SurfaceType::cz, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 2, 16}),
                       this->get_surface_data(SurfaceType::cz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({5}), this->get_surface_data(SurfaceType::pz, 3));
}

TEST_F(RingShapeTest, rotated)
{
    RingShape shape(1.0, 4.0, 1.0, 2.0);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
    static const logic_int expected_logic[] = {
        0l, 1l, LOGIC_NOT, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND, 3l, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-4, -2, -4}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4, -1, 4}), this->bbox.upper());

    static const char* expected_surface_names[] = {"ciz", "coz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(4, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}),
                       this->get_surface_data(SurfaceType::cyo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({16}),
                       this->get_surface_data(SurfaceType::cyo, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2}),
                       this->get_surface_data(SurfaceType::py, 3));
}
