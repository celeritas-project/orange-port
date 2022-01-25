//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/test/tstLogicStack.cc
 * \brief LogicEvaluator class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../LogicStack.hh"

#include "celeritas_test.hh"
#include <iomanip>

using celeritas::detail::LogicStack;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(LogicStackTest, logic_stack)
{
    LogicStack s;
    EXPECT_EQ(0, s.size());
    EXPECT_TRUE(s.empty());

#ifdef REQUIRE_ON
    EXPECT_THROW(s.pop(), assertion);
    EXPECT_THROW(s.top(), assertion);
    EXPECT_THROW(s.apply_not(), assertion);
#endif

    // Add a value [F]
    s.push(false);
    EXPECT_EQ(1, s.size());
    EXPECT_FALSE(s.top());

    for (bool val : {false, true, true, false, false, true, false})
    {
        s.push(val);
        EXPECT_EQ(val, s.top());
    }
    EXPECT_EQ(8, s.size());

    // Reverse order
    int ctr = s.size();
    for (bool val : {false, true, false, false, true, true, false, false})
    {
        EXPECT_EQ(val, s.pop()) << "Failed while popping element " << ctr;
        --ctr;
    }
    EXPECT_EQ(0, ctr);
    EXPECT_EQ(0, s.size());
    EXPECT_TRUE(s.empty());
}

//---------------------------------------------------------------------------//
TEST(LogicStackTest, operators)
{
    LogicStack s;
    s.push(false);
    EXPECT_FALSE(s.top());
    s.apply_not();
    EXPECT_TRUE(s.top());

    s.push(false);
    EXPECT_EQ(2, s.size());
    s.apply_or();
    EXPECT_EQ(true, s.top());
    EXPECT_EQ(1, s.size());

    s.push(true);
    s.apply_and();
    EXPECT_EQ(true, s.top());
    EXPECT_EQ(1, s.size());

    s.push(true);
    s.push(true);
    s.push(false);
    s.apply_and();
    EXPECT_FALSE(s.pop());
    EXPECT_TRUE(s.pop());
    EXPECT_EQ(1, s.size());

    s.push(false);
    s.push(false);
    s.apply_or();
    EXPECT_FALSE(s.top());
    EXPECT_EQ(2, s.size());

    s.apply_or();
    EXPECT_EQ(true, s.top());
    EXPECT_EQ(1, s.size());

    s.push(false);
    s.apply_and();
    EXPECT_FALSE(s.top());
    EXPECT_EQ(1, s.size());
}
