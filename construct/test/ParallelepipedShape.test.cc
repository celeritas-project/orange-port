//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstParallelepipedShape.cc
 * \brief Tests for class ParallelepipedShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../ParallelepipedShape.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::ParallelepipedShape;

using geometria::rotation_matrix;

using Axis::x;
using Axis::y;
using Axis::z;
using constants::pi;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class ParallelepipedShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ParallelepipedShapeTest, cuboid)
{
    const Real3         axis_length(6, 5, 7);
    ParallelepipedShape shape(axis_length, 0, 0, 0);

    EXPECT_EQ("parallelepiped", shape.type());
    EXPECT_SOFT_EQ(210, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);
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

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({6, 5, 7}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({6}), this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({5}), this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({7}), this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ParallelepipedShapeTest, ppiped_psi_angle)
{
    const Real3 axis_length(6, 5, 7);

    ParallelepipedShape shape(axis_length, 55 * pi / 180.0, 0, 0);

    EXPECT_SOFT_EQ(120.45105163371969, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);
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

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({10.095760221445, 2.8678821817552, 7}),
                       this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.57357643635105, -0.81915204428899, 0, 0}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.57357643635105, -0.81915204428899, 0, 3.4414586181063}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({2.8678821817552}),
                       this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({7}), this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ParallelepipedShapeTest, ppiped_theta_angle)
{
    const Real3 axis_length(6, 5, 7);

    ParallelepipedShape shape(axis_length, 0, 55 * pi / 180.0, 0);

    this->build(shape);
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

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({11.734064310023, 5, 4.0150350544573}),
                       this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.57357643635105, 0, -0.81915204428899, 0}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.57357643635105, 0, -0.81915204428899, 3.4414586181063}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({5}), this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({4.0150350544573}),
                       this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ParallelepipedShapeTest, ppiped_phi_angle)
{
    const Real3 axis_length(6, 5, 7);

    ParallelepipedShape shape(axis_length, 0, 0, 55 * pi / 180.0);

    this->build(shape);
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

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({6, 5, 7}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({6}), this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({5}), this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({7}), this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ParallelepipedShapeTest, ppiped_psi_phi_angle)
{
    const Real3 axis_length(6, 5, 7);

    ParallelepipedShape shape(axis_length, 55 * pi / 180.0, 0, 15 * pi / 180.0);

    this->build(shape);
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

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({10.095760221445, 2.8678821817552, 7}),
                       this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.57357643635105, -0.81915204428899, 0, 0}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.57357643635105, -0.81915204428899, 0, 3.4414586181063}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({2.8678821817552}),
                       this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({7}), this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ParallelepipedShapeTest, ppiped_psi_theta_angle)
{
    const Real3 axis_length(6, 5, 7);

    ParallelepipedShape shape(axis_length, 55 * pi / 180.0, 15 * pi / 180.0, 0);

    this->build(shape);
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

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(
        Real3({11.907493537163, 2.8678821817552, 6.7614807840235}),
        this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.56692007005973, -0.8096457680379, -0.1519057749455, 0}),
        this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.56692007005973,
                               -0.8096457680379,
                               -0.1519057749455,
                               3.4015204203584}),
                       this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({2.8678821817552}),
                       this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({6.7614807840235}),
                       this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ParallelepipedShapeTest, ppiped_theta_phi_angle)
{
    const Real3 axis_length(6, 5, 7);

    ParallelepipedShape shape(axis_length, 0, 12 * pi / 180.0, 15 * pi / 180.0);

    this->build(shape);
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

    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(
        Real3({7.4057909022381, 5.3766805369817, 6.8470332051366}),
        this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.97956688411294, 0, -0.201118670315, 0}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.97956688411294, 0, -0.201118670315, 5.8774013046776}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0, 0.99849017341965, -0.054930625195818, 0}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.99849017341965, -0.054930625195818, 4.9924508670982}),
        this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({6.8470332051366}),
                       this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ParallelepipedShapeTest, translated)
{
    const Real3 axis_length(6, 5, 7);

    ParallelepipedShape shape(axis_length, 0, 12 * pi / 180.0, 15 * pi / 180.0);

    this->build(shape, Transform{{1, 2, 3}});
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

    EXPECT_VEC_SOFT_EQ(Real3({1, 2, 3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(
        Real3({8.4057909022381, 7.3766805369817, 9.8470332051366}),
        this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.97956688411294, 0, -0.201118670315, 0.37621087316795}),
        this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.97956688411294, 0, -0.201118670315, 6.2536121778456}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.99849017341965, -0.054930625195818, 1.8321884712518}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.99849017341965, -0.054930625195818, 6.8246393383501}),
        this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::pz, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({9.8470332051366}),
                       this->get_surface_data(SurfaceType::pz, 5));
}

TEST_F(ParallelepipedShapeTest, rotated)
{
    const Real3 axis_length(6, 5, 7);

    ParallelepipedShape shape(axis_length, 0, 12 * pi / 180.0, 15 * pi / 180.0);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});

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
                                               LOGIC_AND,
                                               5l,
                                               LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({0, -6.8470332051366, 0}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({7.4057909022381, 0, 5.3766805369817}),
                       this->bbox.upper());

    static const char* expected_surface_names[]
        = {"p0", "p1", "p2", "p3", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(6, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.97956688411294, 0.201118670315, 0, 0}),
                       this->get_surface_data(SurfaceType::p, 0));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.97956688411294, 0.201118670315, 0, 5.8774013046776}),
        this->get_surface_data(SurfaceType::p, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0, 0.054930625195818, 0.99849017341965, 0}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0, 0.054930625195818, 0.99849017341965, 4.9924508670982}),
        this->get_surface_data(SurfaceType::p, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({-6.8470332051366}),
                       this->get_surface_data(SurfaceType::py, 5));
}
