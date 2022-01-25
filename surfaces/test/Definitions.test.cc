//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/test/tstDefinitions.cc
 * \brief Tests for class Definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Definitions.hh"
#include "../../Definitions.hh"

#include <cmath>
#include <limits>
#include "celeritas_test.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(ConstantsTest, all)
{
    EXPECT_TRUE(std::isinf(no_intersection()));
}

TEST(EnumsTest, sense)
{
    EXPECT_EQ(Sense::neg, to_sense(false));
    EXPECT_EQ(Sense::pos, to_sense(true));
    EXPECT_EQ(Sense::neg, flip_sense(Sense::pos));
    EXPECT_EQ(Sense::pos, flip_sense(Sense::neg));
}

TEST(EnumsTest, signed_sense)
{
    using Limits = std::numeric_limits<real_type>;

    EXPECT_EQ(SignedSense::outside, real_to_sense(Limits::quiet_NaN()));
    EXPECT_EQ(SignedSense::outside, real_to_sense(-Limits::quiet_NaN()));
    EXPECT_EQ(SignedSense::outside, real_to_sense(Limits::infinity()));
    EXPECT_EQ(SignedSense::outside, real_to_sense(1.0));
    EXPECT_EQ(SignedSense::on, real_to_sense(0.0));
    EXPECT_EQ(SignedSense::on, real_to_sense(-0.0));
    EXPECT_EQ(SignedSense::inside, real_to_sense(-1.0));
    EXPECT_EQ(SignedSense::inside, real_to_sense(-Limits::infinity()));

    EXPECT_EQ(Sense::inside, to_sense(SignedSense::inside));   // < 0
    EXPECT_EQ(Sense::pos, to_sense(SignedSense::on));          // >= 0
    EXPECT_EQ(Sense::outside, to_sense(SignedSense::outside)); // >= 0
}

TEST(EnumsTest, zorder)
{
    EXPECT_FALSE(is_exterior_zorder(zorder_int(ZOrder::media)));
    EXPECT_FALSE(is_exterior_zorder(zorder_int(ZOrder::hole)));
    EXPECT_TRUE(is_exterior_zorder(zorder_int(ZOrder::implicit_exterior)));
    EXPECT_TRUE(is_exterior_zorder(zorder_int(ZOrder::exterior)));
}
