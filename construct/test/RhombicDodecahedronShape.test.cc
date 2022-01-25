//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstRhombicDodecahedronShape.cc
 * \brief Tests for class RhombicDodecahedronShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../RhombicDodecahedronShape.hh"

#include "celeritas_test.hh"
#include "orange/TransformUtils.hh"
#include "ShapeTest.hh"

using celeritas::RhombicDodecahedronShape;

using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class RhombicDodecahedronShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(RhombicDodecahedronShapeTest, origin)
{
    const real_type          inscribing_radius = 3;
    RhombicDodecahedronShape shape(inscribing_radius);

    EXPECT_EQ("rhombic_dodecahedron", shape.type());
    EXPECT_SOFT_EQ(152.73506473629428, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(inscribing_radius, shape.inradius());

    this->build(shape);
    EXPECT_VEC_EQ(CSGCell::from_string("0 1 ~ & 2 & 3 ~ & 4 & 5 ~ & 6 ~ & 7 "
                                       "& 8 ~ & 9 & 10 & 11 ~ & ")
                      .logic(),
                  this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-3, -3, -4.24264068711929}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({3, 3, 4.24264068711929}), this->bbox.upper());

    static const char* const expected_surface_names[] = {
        "mx", "px", "my", "py", "ma", "pa", "mb", "pb", "mc", "pc", "md", "pd"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(12, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-3}),
                       this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-3}),
                       this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.5, 0.70710678118655, -3}),
                       this->get_surface_data(SurfaceType::p, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.5, 0.70710678118655, 3}),
                       this->get_surface_data(SurfaceType::p, 5));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.5, -0.70710678118655, 3}),
                       this->get_surface_data(SurfaceType::p, 6));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.5, -0.70710678118655, -3}),
                       this->get_surface_data(SurfaceType::p, 7));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.5, -0.70710678118655, 3}),
                       this->get_surface_data(SurfaceType::p, 8));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.5, -0.70710678118655, -3}),
                       this->get_surface_data(SurfaceType::p, 9));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.5, 0.70710678118655, -3}),
                       this->get_surface_data(SurfaceType::p, 10));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.5, 0.70710678118655, 3}),
                       this->get_surface_data(SurfaceType::p, 11));
}

TEST_F(RhombicDodecahedronShapeTest, translated)
{
    const real_type          inscribing_radius = 1;
    RhombicDodecahedronShape shape(inscribing_radius);

    this->build(shape, Transform{{1, 2, 3}});

    EXPECT_VEC_SOFT_EQ(Real3({0, 1, 1.5857864376269}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 3, 4.4142135623731}), this->bbox.upper());

    ASSERT_EQ(12, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({0}), this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({2}), this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({1}), this->get_surface_data(SurfaceType::py, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({3}), this->get_surface_data(SurfaceType::py, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.5, 0.70710678118655, 2.6213203435596}),
                       this->get_surface_data(SurfaceType::p, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.5, 0.70710678118655, 4.6213203435596}),
                       this->get_surface_data(SurfaceType::p, 5));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.5, -0.70710678118655, -1.6213203435596}),
                       this->get_surface_data(SurfaceType::p, 6));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.5, -0.70710678118655, -3.6213203435596}),
                       this->get_surface_data(SurfaceType::p, 7));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.5, -0.70710678118655, 0.37867965644036}),
                       this->get_surface_data(SurfaceType::p, 8));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.5, -0.70710678118655, -1.6213203435596}),
                       this->get_surface_data(SurfaceType::p, 9));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.5, 0.70710678118655, 0.62132034355964}),
                       this->get_surface_data(SurfaceType::p, 10));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.5, 0.70710678118655, 2.6213203435596}),
                       this->get_surface_data(SurfaceType::p, 11));
}

TEST_F(RhombicDodecahedronShapeTest, rotated)
{
    const real_type          inscribing_radius = 9;
    RhombicDodecahedronShape shape(inscribing_radius);

    // Rotate 90 degrees about X
    this->build(shape, Transform{rotation_matrix(Axis::x, .25)});
    EXPECT_VEC_EQ(CSGCell::from_string("0 1 ~ & 2 & 3 ~ & 4 & 5 ~ & 6 ~ & 7 "
                                       "& 8 ~ & 9 & 10 & 11 ~ & ")
                      .logic(),
                  this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-9, -12.727922061358, -9}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({9, 12.727922061358, 9}), this->bbox.upper());

    ASSERT_EQ(12, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({-9}),
                       this->get_surface_data(SurfaceType::px, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({9}), this->get_surface_data(SurfaceType::px, 1));
    EXPECT_VEC_SOFT_EQ(VecDbl({-9}),
                       this->get_surface_data(SurfaceType::pz, 2));
    EXPECT_VEC_SOFT_EQ(VecDbl({9}), this->get_surface_data(SurfaceType::pz, 3));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.70710678118655, 0.5, -9}),
                       this->get_surface_data(SurfaceType::p, 4));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.70710678118655, 0.5, 9}),
                       this->get_surface_data(SurfaceType::p, 5));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.70710678118655, -0.5, 9}),
                       this->get_surface_data(SurfaceType::p, 6));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.70710678118655, -0.5, -9}),
                       this->get_surface_data(SurfaceType::p, 7));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.70710678118655, 0.5, 9}),
                       this->get_surface_data(SurfaceType::p, 8));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, 0.70710678118655, 0.5, -9}),
                       this->get_surface_data(SurfaceType::p, 9));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.70710678118655, -0.5, -9}),
                       this->get_surface_data(SurfaceType::p, 10));
    EXPECT_VEC_SOFT_EQ(VecDbl({0.5, -0.70710678118655, -0.5, 9}),
                       this->get_surface_data(SurfaceType::p, 11));
}
