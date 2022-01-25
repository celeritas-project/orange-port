//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstWedgeShape.cc
 * \brief Tests for class WedgeShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../WedgeShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::WedgeShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class WedgeShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(WedgeShapeTest, degen_triangle)
{
    WedgeShape shape(1.0, 0.0, 1.0, 1.0);

    EXPECT_EQ("wedge", shape.type());
    EXPECT_SOFT_EQ(0.5, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[] = {0l,
                                               1l,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1, 1}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"my", "p0", "p1", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.70710678118655, 0.70710678118655, 0, 0.70710678118655}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 4));
}

TEST_F(WedgeShapeTest, origin)
{
    WedgeShape shape(4.5, 1.2, 3.2, 6.2);

    EXPECT_EQ("wedge", shape.type());
    EXPECT_SOFT_EQ(44.640000000000008, shape.volume());

    this->build(shape);
    static const logic_int expected_logic[] = {0l,
                                               1l,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4.5, 3.2, 6.2}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"my", "p0", "p1", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.93632917756904, -0.35112344158839, 0, 0}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.69614583607375, 0.71790039345105, 0, 3.1326562623319}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({6.2}),
                       this->get_surface_data(SurfaceType::pz, 4));
}

TEST_F(WedgeShapeTest, translated)
{
    WedgeShape shape(4.5, 1.2, 3.2, 6.2);

    this->build(shape, Transform{{1, 2, 3}});
    static const logic_int expected_logic[] = {0l,
                                               1l,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_NOT,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({1, 2, 3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({5.5, 5.2, 9.2}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"my", "p0", "p1", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::py, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.93632917756904, -0.35112344158839, 0, 0.23408229439226}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.69614583607375, 0.71790039345105, 0, 5.2646028853077}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({9.2}),
                       this->get_surface_data(SurfaceType::pz, 4));
}

TEST_F(WedgeShapeTest, rotated)
{
    WedgeShape shape(4.5, 1.2, 3.2, 6.2);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
    static const logic_int expected_logic[] = {0l,
                                               1l,
                                               LOGIC_AND,
                                               2l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               3l,
                                               LOGIC_NOT,
                                               LOGIC_AND,
                                               4l,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, -6.2, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4.5, 0, 3.2}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"my", "p0", "p1", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.93632917756904, 0, -0.35112344158839, 0}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.69614583607375, 0, 0.71790039345105, 3.1326562623319}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-6.2}),
                       this->get_surface_data(SurfaceType::py, 4));
}
