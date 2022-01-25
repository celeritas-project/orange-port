//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstChordShape.cc
 * \brief Tests for class ChordShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../ChordShape.hh"

#include "celeritas_test.hh"
#include "ShapeTest.hh"

using celeritas::ChordShape;
using celeritas::neg;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class ChordShapeTest : public ::celeritas::ShapeTest
{
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ChordShapeTest, accessors)
{
    ChordShape shape{Axis::y, neg, 10.0};

    // TODO Test shape-specific accessors, and check these ones:
    EXPECT_EQ("chord", shape.type());
    EXPECT_SOFT_EQ(0.0, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());
}

TEST_F(ChordShapeTest, untranformed)
{
    ChordShape shape{Axis::y, neg, 10.0};

    this->build(shape);
    // this->print_expected();

    EXPECT_VEC_EQ(CSGCell::from_string("0 ~ ").logic(), this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({inf, 10, inf}), this->bbox.upper());

    static const char* const expected_surface_names[] = {"py"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({10}),
                       this->get_surface_data(SurfaceType::py, 0));
}
