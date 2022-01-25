//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstGeneralQuadricShape.cc
 * \brief Tests for class GeneralQuadricShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../GeneralQuadricShape.hh"

#include <limits>

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "orange/surfaces/GeneralQuadric.hh"
#include "ShapeTest.hh"

using celeritas::GeneralQuadricShape;

using celeritas::GeneralQuadric;
using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class GeneralQuadricShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(GeneralQuadricShapeTest, accessors)
{
    GeneralQuadricShape shape(1, 1, 1, 0, 0, 0, 0, 0, 0, 0);

    EXPECT_EQ("general_quadric", shape.type());
    EXPECT_SOFT_EQ(0.0, shape.volume());
    EXPECT_FALSE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());
}

//---------------------------------------------------------------------------//
/*!
 * Quadric - origin centered unit sphere
 */
TEST_F(GeneralQuadricShapeTest, sphere)
{
    GeneralQuadricShape pp(1, 1, 1, 0, 0, 0, 0, 0, 0, 0);

    // Build shapes
    this->build(pp);
    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"gq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 1, 1, 0, 0, 0, 0}),
                       this->get_surface_data(SurfaceType::sq, 0));
}

//---------------------------------------------------------------------------//
/*!
 * Quadric - origin centered ellipsoid
 */
TEST_F(GeneralQuadricShapeTest, ellipsoid)
{
    GeneralQuadricShape pp(1, 10, 5, 0, 0, 0, 0, 0, 0, 0);

    // Build shapes
    this->build(pp);

    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"gq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 10, 5, 0, 0, 0, 0}),
                       this->get_surface_data(SurfaceType::sq, 0));
}

//---------------------------------------------------------------------------//
/*!
 * Quadric - 85 degree parallelepiped X planes
 */
TEST_F(GeneralQuadricShapeTest, ppiped_xp_85)
{
    GeneralQuadricShape pp(
        0, 0.00759612, 0.984865, 0, -0.172987, -0, -0, -0.00397227, 0.0452305, 0);

    // Build shapes
    this->build(pp);

    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"gq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0,
                               0.00759612,
                               0.984865,
                               0,
                               -0.172987,
                               0,
                               0,
                               -0.00397227,
                               0.0452305,
                               0}),
                       this->get_surface_data(SurfaceType::gq, 0));
}

//---------------------------------------------------------------------------//
/*!
 * Quadric - yplane
 */
TEST_F(GeneralQuadricShapeTest, yplane)
{
    GeneralQuadricShape pp(0, 0, 0, 0, 0, 0, 0, 1, 0, 0);

    // Build shapes
    this->build(pp);

    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, 0, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"gq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::py, 0));
}

//---------------------------------------------------------------------------//
/*!
 * Quadric - xplane
 */
TEST_F(GeneralQuadricShapeTest, xplane)
{
    GeneralQuadricShape pp(0, 0, 0, 0, 0, 0, -1, 0, 0, 2);

    // Build shapes
    this->build(pp);

    static const logic_int expected_logic[] = {0l};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({2, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"gq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::px, 0));
}

//---------------------------------------------------------------------------//
/*!
 * Quadric - zplane
 */
TEST_F(GeneralQuadricShapeTest, zplane)
{
    GeneralQuadricShape pp(0, 0, 0, 0, 0, 0, 0, 0, 1, 3);

    // Build shapes
    this->build(pp);

    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, inf, -3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"gq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-3}),
                       this->get_surface_data(SurfaceType::pz, 0));
}

//---------------------------------------------------------------------------//
/*!
 * Quadric - zplane
 */
TEST_F(GeneralQuadricShapeTest, zplane_xydiag)
{
    GeneralQuadricShape pp(0, 0, 0, 0, 0, 0, 1, 1, 0, 0);

    // Unknown
    ASSERT_EQ(0, pp.volume());
    ASSERT_EQ(0, pp.inradius());

    // Build shapes
    this->build(pp);

    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), this->bbox.upper());

    static const char* expected_surface_names[] = {"gq"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.70710678118655, 0.70710678118655, 0, 0}),
                       this->get_surface_data(SurfaceType::p, 0));
}
