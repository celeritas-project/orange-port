//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstTriangularPrismShape.cc
 * \brief Tests for class TriangularPrismShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../TriangularPrismShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::TriangularPrismShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class TriangularPrismShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(TriangularPrismShapeTest, origin)
{
    TriangularPrismShape shape(1, 2, 3);

    EXPECT_EQ("triangular_prism", shape.type());
    EXPECT_SOFT_EQ(3.0, shape.volume());
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

    EXPECT_VEC_SOFT_EQ(Real3({-0.5, -1, -1.5}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({0.5, 1, 1.5}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"mx", "my", "p0", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-0.5}),
                       this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.89442719099992, 0.44721359549996, 0, 0}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1.5}),
                       this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({1.5}),
                       this->get_surface_data(SurfaceType::pz, 4));
}

TEST_F(TriangularPrismShapeTest, translated)
{
    TriangularPrismShape shape(1, 2, 3);

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

    EXPECT_VEC_SOFT_EQ(Real3({0.5, 1, 1.5}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1.5, 3, 4.5}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"mx", "my", "p0", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5}),
                       this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.89442719099992, 0.44721359549996, 0, 1.7888543819998}),
        this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({1.5}),
                       this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({4.5}),
                       this->get_surface_data(SurfaceType::pz, 4));
}

TEST_F(TriangularPrismShapeTest, rotated)
{
    TriangularPrismShape shape(1, 2, 3);

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

    EXPECT_VEC_SOFT_EQ(Real3({-0.5, -1.5, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({0.5, 1.5, 1}), this->bbox.upper());

    static const char* expected_surface_names[]
        = {"mx", "my", "p0", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(5, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-0.5}),
                       this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.89442719099992, 0, 0.44721359549996, 0}),
                       this->get_surface_data(SurfaceType::p, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({1.5}),
                       this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1.5}),
                       this->get_surface_data(SurfaceType::py, 4));
}
