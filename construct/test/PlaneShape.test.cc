//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstPlaneShape.cc
 * \brief Tests for class PlaneShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../PlaneShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "orange/surfaces/Plane.hh"
#include "ShapeTest.hh"

using celeritas::PlaneShape;

using celeritas::Plane;
using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class PlaneShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(PlaneShapeTest, accessors)
{
    PlaneShape shape(Real3{1, 2, 0}, Real3{0, 1, 0});

    EXPECT_VEC_SOFT_EQ(Real3(1, 2, 0), shape.normal());
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), shape.point());
    EXPECT_EQ("plane", shape.type());
    EXPECT_SOFT_EQ(0.0, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());
}

TEST_F(PlaneShapeTest, from_displacement)
{
    PlaneShape shape = PlaneShape::from_displacement(Real3{10, 0, 0}, 5.0);
    EXPECT_VEC_SOFT_EQ(Real3(10, 0, 0), shape.normal());
    EXPECT_VEC_SOFT_EQ(Real3(0.5, 0, 0), shape.point());
    this->build(shape);

    EXPECT_VEC_EQ(CSGCell::from_string("0 ~ ").logic(), this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({0.5, inf, inf}), this->bbox.upper());

    static const char* const expected_surface_names[] = {"p0"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5}),
                       this->get_surface_data(SurfaceType::px, 0));
}

TEST_F(PlaneShapeTest, untransformed)
{
    PlaneShape shape(Real3{1, 2, 0}, Real3{0, 1, 0});

    this->build(shape);

    EXPECT_VEC_EQ(CSGCell::from_string("0 ~ ").logic(), this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), this->bbox.upper());

    static const char* const expected_surface_names[] = {"p0"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(
        VecDbl({0.44721359549996, 0.89442719099992, 0, 0.89442719099992}),
        this->get_surface_data(SurfaceType::p, 0));
}

TEST_F(PlaneShapeTest, translated)
{
    PlaneShape shape(Real3{1, 0, 0}, Real3{10, 0, 0});

    this->build(shape, Transform{{1, 2, 3}});

    EXPECT_VEC_EQ(CSGCell::from_string("0 ~ ").logic(), this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({11, inf, inf}), this->bbox.upper());

    static const char* const expected_surface_names[] = {"p0"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({11}),
                       this->get_surface_data(SurfaceType::px, 0));
}

TEST_F(PlaneShapeTest, rotated)
{
    PlaneShape shape(Real3{0, 1, 0}, Real3{3, 0, 0});

    {
        // Rotate 45 degrees
        this->build(shape, Transform{rotation_matrix(Axis::x, .125)});
        EXPECT_VEC_EQ(CSGCell::from_string("0 ~ ").logic(), this->cell.logic());

        EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), this->bbox.upper());

        static const char* const expected_surface_names[] = {"p0"};
        EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

        ASSERT_EQ(1, this->surfaces.size());
        EXPECT_VEC_SOFT_EQ(VecDbl({0, 0.70710678118655, 0.70710678118655, 0}),
                           this->get_surface_data(SurfaceType::p, 0));
    }
    {
        // Rotate 90 degrees about X
        this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
        EXPECT_VEC_EQ(CSGCell::from_string("0 ~ ").logic(), this->cell.logic());

        EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({inf, inf, 0}), this->bbox.upper());

        static const char* const expected_surface_names[] = {"p0"};
        EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

        ASSERT_EQ(1, this->surfaces.size());
        EXPECT_VEC_SOFT_EQ(VecDbl({0}),
                           this->get_surface_data(SurfaceType::pz, 0));
    }
    {
        // Rotate 180 degrees about X: note that the CSG cell logic is negated
        this->build(shape, Transform{rotation_matrix(Axis::x, .5)});
        EXPECT_VEC_EQ(CSGCell::from_string("0 ").logic(), this->cell.logic());

        EXPECT_VEC_SOFT_EQ(Real3({-inf, 0, -inf}), this->bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), this->bbox.upper());

        static const char* const expected_surface_names[] = {"p0"};
        EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

        ASSERT_EQ(1, this->surfaces.size());
        EXPECT_VEC_SOFT_EQ(VecDbl({0}),
                           this->get_surface_data(SurfaceType::py, 0));
    }
}
