//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstCuboidShape.cc
 * \brief Tests for class CuboidShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../CuboidShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "orange/surfaces/PlaneAligned.hh"
#include "orange/surfaces/Plane.hh"
#include "ShapeTest.hh"

using celeritas::CuboidShape;
using celeritas::Plane;
using celeritas::PlaneX;
using celeritas::PlaneY;
using celeritas::PlaneZ;
using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class CuboidShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CuboidShapeTest, origin)
{
    Real3       lo{-3, -2, 0};
    Real3       hi{3, 1, 4};
    CuboidShape shape(lo, hi);

    EXPECT_VEC_SOFT_EQ(lo, shape.lower());
    EXPECT_VEC_SOFT_EQ(hi, shape.upper());
    EXPECT_EQ("cuboid", shape.type());
    EXPECT_SOFT_EQ(6.0 * 3.0 * 4.0, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(1.5, shape.inradius());

    this->build(shape);

    static const logic_int expected_logic[] = {0ul,
                                               1ul,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2ul,
                                               LOGIC_AND,
                                               3ul,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               4ul,
                                               LOGIC_AND,
                                               5ul,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    static const real_type expected_bbox_lower[] = {-3, -2, 0};
    static const real_type expected_bbox_upper[] = {3, 1, 4};
    static const char*     expected_surface_names[]
        = {"mx", "px", "my", "py", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    EXPECT_VEC_SOFT_EQ(lo, this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(hi, this->bbox.upper());

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_SOFT_EQ(-3.0, this->get_surface<PlaneX>(0).position());
    EXPECT_SOFT_EQ(3.0, this->get_surface<PlaneX>(1).position());
    EXPECT_SOFT_EQ(-2.0, this->get_surface<PlaneY>(2).position());
    EXPECT_SOFT_EQ(1.0, this->get_surface<PlaneY>(3).position());
    EXPECT_SOFT_EQ(0.0, this->get_surface<PlaneZ>(4).position());
    EXPECT_SOFT_EQ(4.0, this->get_surface<PlaneZ>(5).position());
}

//---------------------------------------------------------------------------//

TEST_F(CuboidShapeTest, translation)
{
    CuboidShape shape({0, 0, 0}, {1, 1, 1});

    this->build(shape, Transform{Real3{10, 0, 0}});
    EXPECT_VEC_SOFT_EQ(Real3({10, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({11, 1, 1}), this->bbox.upper());
}

//---------------------------------------------------------------------------//
//! Rotate off the main axis. Plane -y sense should be flipped.

TEST_F(CuboidShapeTest, simple_rotation)
{
    Real3       lo{0, 0, 0};
    Real3       hi{4, 3, 4};
    CuboidShape shape(lo, hi);

    this->build(shape, Transform{rotation_matrix(Axis::z, 1. / 6)});

    // Check surface data
    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.86602540378444, 0, 0}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.86602540378444, 0, 4}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, 0}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, -3}),
                       this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({4}), this->get_surface_data(SurfaceType::pz, 5));
}

//---------------------------------------------------------------------------//
//! Rotate off the main axis. Plane -y sense should be flipped.
TEST_F(CuboidShapeTest, off_axis_rotation)
{
    CuboidShape shape
        = CuboidShape::from_bounds({1.0, 5.0}, {1.0, 4.0}, {0.0, 4.0});

    this->build(shape, Transform{rotation_matrix(Axis::z, 1. / 6)});

    // Logic for surfaces 2 and 3 should be negated
    static const logic_int expected_logic[] = {0ul,
                                               1ul,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               2ul,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               3ul,
                                               LOGIC_AND,
                                               4ul,
                                               LOGIC_AND,
                                               5ul,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    // Check surface data
    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.86602540378444, 0, 1}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.86602540378444, 0, 5}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, -1}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.86602540378444, -0.5, 0, -4}),
                       this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({4}), this->get_surface_data(SurfaceType::pz, 5));
}

//---------------------------------------------------------------------------//
//! Rotate about z, y, and x so that all senses should be flipped.
TEST_F(CuboidShapeTest, upside_down)
{
    Real3       lo{-3, -2, -1};
    Real3       hi{4, 5, 6};
    CuboidShape shape(lo, hi);

    this->build(
        shape,
        Transform{rotation_matrix(Axis::z, 0.5, Axis::y, 0.5, Axis::x, 0.5)});

    // this->print_expected();

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
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-3, -2, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4, 5, 6}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"mx", "px", "my", "py", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-3}),
                       this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({4}), this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-2}),
                       this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({5}), this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({6}), this->get_surface_data(SurfaceType::pz, 5));
}
