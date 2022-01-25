//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstBBoxCalculator.cc
 * \brief Tests for class BBoxCalculator
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../BBoxCalculator.hh"

#include <limits>
#include "celeritas_test.hh"
#include "../SlabShape.hh"
#include "../CuboidShape.hh"
#include "../PlacedShape.hh"

using namespace celeritas;

static const real_type inf = std::numeric_limits<real_type>::infinity();

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(BBoxCalculatorTest, all)
{
    BBoxCalculator calc_bbox;
    {
        auto bbox = calc_bbox(SlabShape(Axis::z, -3, 3));

        EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -3}), bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({inf, inf, 3}), bbox.upper());
    }

    {
        auto bbox = calc_bbox(CuboidShape({0, 0, 0}, {1, 2, 3}),
                              Transform({1, 2, 3}));

        EXPECT_VEC_SOFT_EQ(Real3({1, 2, 3}), bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({2, 4, 6}), bbox.upper());
    }

    {
        PlacedShape::Params params;
        params.shape     = std::make_shared<SlabShape>(Axis::x, 0, 2);
        params.transform = Transform({1, 0, 0});
        params.md        = ORANGE_MD_FROM_SOURCE("test-shape");
        PlacedShape ps(std::move(params));

        auto bbox = calc_bbox(ps);

        EXPECT_VEC_SOFT_EQ(Real3({1, -inf, -inf}), bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({3, inf, inf}), bbox.upper());
    }

    // Surface deduplication should not be taking place with the current
    // implementation; however, we might be able to speed up bbox calculation
    // for lots of shapes if we reuse the surfaces container and surface
    // inserter.
    {
        auto bbox
            = calc_bbox(SlabShape(Axis::x, 0, 2), Transform({1e-9, 0, 0}));

        EXPECT_VEC_SOFT_EQ(Real3({1e-9, -inf, -inf}), bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({2 + 1e-9, inf, inf}), bbox.upper());
    }
}
