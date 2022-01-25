//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstConeShape.cc
 * \brief Tests for class ConeShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../ConeShape.hh"

#include <cmath>

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/Definitions.hh"
#include "orange/TransformUtils.hh"
#include "orange/surfaces/ConeAligned.hh"
#include "ShapeTest.hh"

using celeritas::ConeShape;

using celeritas::ConeAligned;
using constants::pi;
using geometria::Real3;
using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class ConeShapeTest : public ::celeritas::ShapeTest
{
  protected:
    //-----------------------------------------------------------------------//
    void test_accessors(const ConeShape& shape,
                        real_type        hr,
                        real_type        lr,
                        real_type        h,
                        real_type        l)
    {
        // Expected inscribing radius
        real_type mr = std::max(hr, lr); // max radius indicates max x, y
        real_type hh = 0.5 * (h - l);    // half height
        // Sqrt of the sum of the squares in distance from the origin
        // to the max x, y, and z
        real_type expected_ir = std::sqrt(2.0 * ipow<2>(mr) + ipow<2>(hh));

        // Expected volume is pi/3 * H * lr^2 + hr^2 + lr*hr
        real_type expected_vol
            = (h - l) * (ipow<2>(lr) + lr * hr + ipow<2>(hr)) * pi / 3;

        // Test shape-specific accessors, and check these ones:
        EXPECT_EQ("cone", shape.type());
        EXPECT_SOFT_EQ(expected_vol, shape.volume());
        EXPECT_TRUE(shape.is_convex());
        EXPECT_SOFT_EQ(expected_ir, shape.inradius());

        EXPECT_EQ(hr, shape.high_radius());
        EXPECT_EQ(lr, shape.low_radius());
        EXPECT_EQ(h, shape.high_extent());
        EXPECT_EQ(l, shape.low_extent());
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ConeShapeTest, origin)
{
    // Radii
    real_type hr = 1.2; // high radius
    real_type lr = 0.0; // low radius

    // Axis intersects
    real_type h = 1.1;  // high
    real_type l = -0.2; // low

    ConeShape shape(Axis::z, {lr, hr}, {l, h});

    this->test_accessors(shape, hr, lr, h, l);

    this->build(shape);

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1.2, -1.2, -0.2}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1.2, 1.2, 1.1}), this->bbox.upper());

    static const char* expected_surface_names[] = {"koz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0, 0, -0.2, 0.85207100591716}),
                       this->get_surface_data(SurfaceType::kz, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-0.2}),
                       this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1.1}),
                       this->get_surface_data(SurfaceType::pz, 2));
}

TEST_F(ConeShapeTest, translated)
{
    // Radii
    real_type hr = 1.2; // high radius
    real_type lr = 1.3; // low radius

    // Axis intersects
    real_type h = 0.7;  // high
    real_type l = -0.7; // low

    ConeShape shape(Axis::z, {lr, hr}, {l, h});

    this->test_accessors(shape, hr, lr, h, l);

    Real3 t = {1.0, 2.0, 3.0};
    this->build(shape, Transform{t});

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_AND, 2l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-0.3, 0.7, 2.3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2.3, 3.3, 3.7}), this->bbox.upper());

    static const char* expected_surface_names[] = {"koz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 2, 20.5, 0.0051020408163265}),
                       this->get_surface_data(SurfaceType::kz, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({2.3}),
                       this->get_surface_data(SurfaceType::pz, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({3.7}),
                       this->get_surface_data(SurfaceType::pz, 2));
}

TEST_F(ConeShapeTest, rotated)
{
    // Radii
    real_type hr = 1.2; // high radius
    real_type lr = 1.3; // low radius

    // Axis intersects
    real_type h = 0.7;  // high
    real_type l = -0.7; // low

    ConeShape shape(Axis::z, {lr, hr}, {l, h});

    this->test_accessors(shape, hr, lr, h, l);

    // Rotate 90 degrees about X
    Transform t = {rotation_matrix(Axis::x, .25)};
    this->build(shape, t);

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_NOT, LOGIC_AND, 2l, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-1.3, -0.7, -1.3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({1.3, 0.7, 1.3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"koz", "mz", "pz"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(3, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-0, -17.5, -0, 0.0051020408163265}),
                       this->get_surface_data(SurfaceType::ky, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.7}),
                       this->get_surface_data(SurfaceType::py, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-0.7}),
                       this->get_surface_data(SurfaceType::py, 2));
}
