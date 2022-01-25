//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstSlabShape.cc
 * \brief Tests for class SlabShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../SlabShape.hh"

#include <limits>

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::SlabShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class SlabShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(SlabShapeTest, origin)
{
    SlabShape shape(Axis::z, -3, 3);

    EXPECT_EQ("slab", shape.type());
    EXPECT_SOFT_EQ(0.0, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(3.0, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {0l, 1l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, inf, 3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(2, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-3}),
                       this->get_surface_data(SurfaceType::pz, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::pz, 1));
}

TEST_F(SlabShapeTest, translated)
{
    SlabShape shape(Axis::x, -3, 5);

    this->build(shape, Transform{{1, 2, 3}});
    static const logic_int expected_logic[] = {0l, 1l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({6, inf, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"mx", "px"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(2, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-2}),
                       this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({6}), this->get_surface_data(SurfaceType::px, 1));
}

TEST_F(SlabShapeTest, rotated)
{
    SlabShape shape(Axis::z, -3, 5);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
    static const logic_int expected_logic[] = {0l, LOGIC_NOT, 1l, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -5, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, 3, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(2, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::py, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-5}),
                       this->get_surface_data(SurfaceType::py, 1));
}
