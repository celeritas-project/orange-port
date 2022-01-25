//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstCSGCell.cc
 * \brief Tests for class CSGCell
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../CSGCell.hh"

#include "celeritas_test.hh"

using celeritas::CSGCell;
using celeritas::SurfaceId;
using logic_int = CSGCell::logic_int;

static constexpr auto LOGIC_TRUE = CSGCell::LOGIC_TRUE;
static constexpr auto LOGIC_OR   = CSGCell::LOGIC_OR;
static constexpr auto LOGIC_AND  = CSGCell::LOGIC_AND;
static constexpr auto LOGIC_NOT  = CSGCell::LOGIC_NOT;

std::vector<int> vec_surfid_to_int(const std::vector<SurfaceId>& inp)
{
    std::vector<int> result;
    for (SurfaceId id : inp)
    {
        result.push_back(id.unchecked_get());
    }
    return result;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(UnitRegionTest, logic_values)
{
    EXPECT_GE(logic_int(LOGIC_TRUE), CSGCell::BEGIN_LOGIC_TOKEN);
    EXPECT_GE(logic_int(LOGIC_OR), CSGCell::BEGIN_LOGIC_TOKEN);
    EXPECT_GE(logic_int(LOGIC_AND), CSGCell::BEGIN_LOGIC_TOKEN);
    EXPECT_GE(logic_int(LOGIC_NOT), CSGCell::BEGIN_LOGIC_TOKEN);
    EXPECT_TRUE(CSGCell::is_operator_token(LOGIC_TRUE));
    EXPECT_TRUE(CSGCell::is_operator_token(LOGIC_OR));
    EXPECT_TRUE(CSGCell::is_operator_token(LOGIC_AND));
    EXPECT_TRUE(CSGCell::is_operator_token(LOGIC_NOT));
    EXPECT_FALSE(CSGCell::is_operator_token(0));
    EXPECT_FALSE(CSGCell::is_operator_token(1));

    EXPECT_EQ('*', to_char(LOGIC_TRUE));
    EXPECT_EQ('|', to_char(LOGIC_OR));
    EXPECT_EQ('&', to_char(LOGIC_AND));
    EXPECT_EQ('~', to_char(LOGIC_NOT));
}

TEST(UnitRegionTest, alpha)
{
    const logic_int alpha_logic[] = {1,
                                     4,
                                     LOGIC_NOT,
                                     LOGIC_AND,
                                     3,
                                     LOGIC_AND,
                                     2,
                                     LOGIC_NOT,
                                     LOGIC_AND,
                                     LOGIC_NOT,
                                     LOGIC_NOT,
                                     8,
                                     LOGIC_NOT,
                                     LOGIC_NOT,
                                     LOGIC_AND,
                                     LOGIC_NOT};

    auto url = CSGCell::from_string("1 4 ~ & 3 & 2 ~ & ~ ~ 8 ~ ~ & ~");
    EXPECT_VEC_EQ(alpha_logic, url.logic());

    EXPECT_EQ(2, url.logic_depth());

    const int expected[] = {1, 2, 3, 4, 8};
    EXPECT_VEC_EQ(expected, vec_surfid_to_int(url.get_local_surfids()));
}

TEST(UnitRegionTest, gamma)
{
    auto url = CSGCell::from_string("8");
    EXPECT_EQ(1, url.logic_depth());

    const int expected[] = {8};
    EXPECT_VEC_EQ(expected, vec_surfid_to_int(url.get_local_surfids()));
}
