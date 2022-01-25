//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstSphereShape.cc
 * \brief Tests for class SphereShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../SphereShape.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/surfaces/CenteredSphere.hh"
#include "orange/surfaces/Sphere.hh"
#include "ShapeTest.hh"

using celeritas::SphereShape;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class SphereShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(SphereShapeTest, origin)
{
    SphereShape shape(2.5);
    EXPECT_SOFT_EQ(2.5, shape.radius());

    EXPECT_EQ("sphere", shape.type());
    EXPECT_SOFT_EQ(std::pow(2.5, 3) * constants::pi * 4.0 / 3.0,
                   shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(2.5, shape.inradius());

    this->build(shape);

    static const logic_int expected_logic[]         = {0, LOGIC_NOT};
    static const real_type expected_bbox_lower[]    = {-2.5, -2.5, -2.5};
    static const real_type expected_bbox_upper[]    = {2.5, 2.5, 2.5};
    static const char*     expected_surface_names[] = {"s"};

    EXPECT_VEC_SOFT_EQ(expected_bbox_lower, this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(expected_bbox_upper, this->bbox.upper());
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_SOFT_EQ(2.5 * 2.5,
                   this->get_surface<celeritas::CenteredSphere>(0).radius_sq());
}

TEST_F(SphereShapeTest, translated)
{
    SphereShape shape(3.0);

    this->build(shape, Transform{{10, 5, 2}});

    static const logic_int expected_logic[]         = {0, LOGIC_NOT};
    static const real_type expected_bbox_lower[]    = {7, 2, -1};
    static const real_type expected_bbox_upper[]    = {13, 8, 5};
    static const char*     expected_surface_names[] = {"s"};

    EXPECT_VEC_EQ(expected_logic, this->cell.logic());
    EXPECT_VEC_SOFT_EQ(expected_bbox_lower, this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(expected_bbox_upper, this->bbox.upper());
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_SOFT_EQ(9.0, this->get_surface<celeritas::Sphere>(0).radius_sq());
    EXPECT_VEC_SOFT_EQ(Real3({10, 5, 2}),
                       this->get_surface<celeritas::Sphere>(0).origin());
}
