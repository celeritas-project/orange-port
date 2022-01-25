//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstECylinderShape.cc
 * \brief Tests for class ECylinderShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../ECylinderShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::ECylinderShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class ECylinderShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ECylinderShapeTest, origin)
{
    ECylinderShape shape({1, 1}, -1, 1);

    EXPECT_EQ("ecylinder", shape.type());
    EXPECT_SOFT_EQ(6.2831853071795862, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(1.7320508075688772, shape.inradius());

    this->build(shape);
    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1, -1, -1}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 1, 1}), this->bbox.upper());

    static const char* expected_surface_names[] = {"sq", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1}),
                       this->get_surface_data(SurfaceType::czo, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1}),
                       this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::pz, 2));
}

TEST_F(ECylinderShapeTest, translated)
{
    ECylinderShape shape({1.1, 2.2}, 0, 2);
    EXPECT_SOFT_EQ(15.205308443374602, shape.volume());
    EXPECT_SOFT_EQ(3.2680269276736387, shape.inradius());

    this->build(shape, Transform{{1, 2, 3}});
    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-0.1, -0.2, 3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2.1, 4.2, 5}), this->bbox.upper());

    static const char* expected_surface_names[] = {"sq", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({4.84, 1.21, 0, -9.68, -4.84, 0, 3.8236}),
                       this->get_surface_data(SurfaceType::sq, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({5}), this->get_surface_data(SurfaceType::pz, 2));
}

TEST_F(ECylinderShapeTest, rotated)
{
    ECylinderShape shape({0.2, 3}, -2, -1);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_NOT, LOGIC_AND, 2l, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-0.2, 1, -3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({0.2, 2, 3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"sq", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({9, 0, 0.04, 0, 0, 0, -0.36}),
                       this->get_surface_data(SurfaceType::sq, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::py, 2));
}
