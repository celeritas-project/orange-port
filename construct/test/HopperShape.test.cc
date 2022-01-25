//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstHopperShape.cc
 * \brief Tests for class HopperShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../HopperShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::HopperShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class HopperShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(HopperShapeTest, big_top)
{
    HopperShape shape(2.0,
                      2.5,
                      4.0, // hi half-widths and z point
                      1.0,
                      1.5,
                      3.0); // lo half-widths and z

    EXPECT_EQ("hopper", shape.type());
    EXPECT_SOFT_EQ(12.318150383367774, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(3.2403703492039302, shape.inradius());

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

    EXPECT_VEC_SOFT_EQ(Real3({-2, -2.5, 3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 2.5, 4}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, 0, -0.70710678118655, -1.4142135623731}),
        this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.70710678118655, -0.70710678118655, -1.0606601717798}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, 0, 0.70710678118655, 1.4142135623731}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.70710678118655, 0.70710678118655, 1.0606601717798}),
        this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({4}), this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(HopperShapeTest, little_top)
{
    HopperShape shape(1.0,
                      1.5,
                      4.0, // hi half-widths and z point
                      2.0,
                      2.5,
                      3.0); // lo half-widths and z

    EXPECT_EQ("hopper", shape.type());
    EXPECT_SOFT_EQ(12.318150383367774, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(3.2403703492039302, shape.inradius());

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

    EXPECT_VEC_SOFT_EQ(Real3({-2, -2.5, 3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 2.5, 4}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, 0, 0.70710678118655, 3.5355339059327}),
        this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.70710678118655, 0.70710678118655, 3.889087296526}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, 0, -0.70710678118655, -3.5355339059327}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.70710678118655, -0.70710678118655, -3.889087296526}),
        this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({4}), this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(HopperShapeTest, translated)
{
    HopperShape shape(1.0,
                      1.5,
                      4.0, // hi half-widths and z point
                      2.0,
                      2.5,
                      3.0); // lo half-widths and z

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
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -0.5, 6}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({3, 4.5, 7}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, 0, 0.70710678118655, 6.3639610306789}),
        this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.70710678118655, 0.70710678118655, 7.4246212024587}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, 0, -0.70710678118655, -4.9497474683058}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.70710678118655, -0.70710678118655, -4.5961940777126}),
        this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({6}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({7}), this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(HopperShapeTest, rotated)
{
    HopperShape shape(1.0,
                      1.5,
                      4.0, // hi half-widths and z point
                      2.0,
                      2.5,
                      3.0); // lo half-widths and z

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
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
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -4, -2.5}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, -3, 2.5}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, -0.70710678118655, 0, 3.5355339059327}),
        this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.70710678118655, -0.70710678118655, -3.889087296526}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, 0.70710678118655, 0, -3.5355339059327}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.70710678118655, 0.70710678118655, -3.889087296526}),
        this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-3}),
                       this->get_surface_data(SurfaceType::py, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({-4}),
                       this->get_surface_data(SurfaceType::py, 5));
}
