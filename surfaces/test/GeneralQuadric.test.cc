//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstGeneralQuadric.cc
 * \brief GeneralQuadric class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../GeneralQuadric.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/Transform.hh"
#include "orange/TransformUtils.hh"

using celeritas::GeneralQuadric;
using celeritas::Real3;
using celeritas::Transform;
using constants::pi;
using geometria::rotation_matrix;

const real_type sqrt_half = std::sqrt(0.5);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(GeneralQuadricTest, construction)
{
    /* From the MCNP manual example:
     * This is a cylinder of radius 1 cm whose axis is in a plane normal to the
     * x–axis at x = 6, displaced 2 cm from the x–axis and rotated 30° about
     * the x–axis off the y–axis toward the z–axis. The sense is positive for
     * points outside the cylinder.
     */
    GeneralQuadric gq({1, .25, .75}, {0, -.866, 0}, {-12, -2, 3.464}, 39);
    EXPECT_VEC_SOFT_EQ(Real3(1, .25, .75), gq.second());
    EXPECT_VEC_SOFT_EQ(Real3(0, -.866, 0), gq.cross());
    EXPECT_VEC_SOFT_EQ(Real3(-12, -2, 3.464), gq.first());
    EXPECT_SOFT_EQ(39, gq.zeroth());

    cout << gq << "\n";

    // Rotate 90 degrees clockwise around +Z,
    Transform t(rotation_matrix(Axis::z, .25), Real3(0, 0, 0));

    auto transformed = gq.transformed(t);

    EXPECT_VEC_SOFT_EQ(Real3(.25, 1, .75), transformed.second());
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, .866), transformed.cross());
    EXPECT_VEC_SOFT_EQ(Real3(2, -12, 3.464), transformed.first());
    EXPECT_SOFT_EQ(39, transformed.zeroth());
}

TEST(PlaneTest, plane_translation)
{
    // GQ for a plane at X=2
    GeneralQuadric gq({0, 0, 0}, {0, 0, 0}, {1, 0, 0}, -2);
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 0), gq.second());
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 0), gq.cross());
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), gq.first());
    EXPECT_SOFT_EQ(-2, gq.zeroth());

    // Move +1 in X; shouldn't modify Z at all
    Transform t(Real3(1, 0, -1));

    auto translated = gq.translated(t);

    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 0), translated.second());
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 0), translated.cross());
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), translated.first());
    EXPECT_SOFT_EQ(-3, translated.zeroth());
}

TEST(PlaneTest, plane_transformation)
{
    // GQ for a plane at X=2
    GeneralQuadric gq({0, 0, 0}, {0, 0, 0}, {1, 0, 0}, -2);

    // Rotate 45 degrees clockwise around +Z, then move +1 in X;
    // Shouldn't modify Z at all
    Transform t(rotation_matrix(Axis::z, .125), Real3(1, 0, -1));

    auto transformed = gq.transformed(t);

    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 0), transformed.second());
    EXPECT_VEC_SOFT_EQ(Real3(0, 0, 0), transformed.cross());
    EXPECT_VEC_SOFT_EQ(Real3(sqrt_half, sqrt_half, 0), transformed.first());
    EXPECT_SOFT_EQ(-2 - sqrt_half, transformed.zeroth());
}
