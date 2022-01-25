//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/test/tstLogicEvaluator.cc
 * \brief LogicEvaluator class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../LogicEvaluator.hh"

#include "celeritas_test.hh"
#include <iomanip>

using make_span;
using celeritas::detail::LogicEvaluator;

using SenseContainer      = celeritas::SenseContainer;
using logic_int           = LogicEvaluator::logic_int;
constexpr auto logic_begin  = LogicEvaluator::logic_begin;
constexpr auto logic_true = LogicEvaluator::logic_true;
constexpr auto logic_or   = LogicEvaluator::logic_or;
constexpr auto logic_and  = LogicEvaluator::logic_and;
constexpr auto logic_not  = LogicEvaluator::logic_not;
constexpr auto logic_end  = LogicEvaluator::logic_end;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(LogicEvaluatorTest, enumeration)
{
    cout << "Begin token: " << std::hex << static_cast<long>(logic_begin)
         << ", end token: " << std::hex << static_cast<long>(logic_end)
         << "\n";

    EXPECT_GE(logic_not, logic_begin);
    EXPECT_GE(logic_and, logic_begin);
    EXPECT_GE(logic_or, logic_begin);

    EXPECT_LT(logic_begin, logic_end);
}

//---------------------------------------------------------------------------//
//! Same logic evaluations as rpn.py

TEST(LogicEvaluatorTest, evaluate)
{
    // Logic for alpha : 1 2 ~ & 3 & 4 ~ & ~ ~ 8 ~ ~ & ~
    // With senses substituted: T F ~ & T & F ~ & T & ~
    const logic_int alpha_logic[] = {1,
                                     2,
                                     logic_not,
                                     logic_and,
                                     3,
                                     logic_and,
                                     4,
                                     logic_not,
                                     logic_and,
                                     logic_not,
                                     logic_not,
                                     8,
                                     logic_not,
                                     logic_not,
                                     logic_and,
                                     logic_not};

    // Logic for beta : 5 1 ~ & 6 & 7 ~ & ~ ~ 8 ~ ~ & ~
    // With senses substituted: T T ~ & F & F ~ & T & ~
    const logic_int beta_logic[] = {5,
                                    1,
                                    logic_not,
                                    logic_and,
                                    6,
                                    logic_and,
                                    7,
                                    logic_not,
                                    logic_and,
                                    logic_not,
                                    logic_not,
                                    8,
                                    logic_not,
                                    logic_not,
                                    logic_and,
                                    logic_not};

    // Logic for gamma : 8 ~ ~ ~ ~
    // With senses substituted: T
    const logic_int gamma_logic[] = {8};

    // 1 2 ~ & 3 & 4 ~ & ~ 5 1 ~ & 6 & 7 ~ & ~ & 8 & 0 ~ &
    const logic_int delta_logic[]
        = {1,         2,         logic_not, logic_and, 3,         logic_and,
           4,         logic_not, logic_and, logic_not, 5,         1,
           logic_not, logic_and, 6,         logic_and, 7,         logic_not,
           logic_and, logic_not, logic_and, 8,         logic_and, 0,
           logic_not, logic_and};

    const logic_int everywhere_logic[] = {logic_true};

    //// CREATE ////

    LogicEvaluator eval_alpha(make_span(alpha_logic));
    LogicEvaluator eval_beta(make_span(beta_logic));
    LogicEvaluator eval_gamma(make_span(gamma_logic));
    LogicEvaluator eval_delta(make_span(delta_logic));
    LogicEvaluator eval_everywhere(make_span(everywhere_logic));

    //// EVALUATE ////

    SenseContainer senses = {0, 1, 0, 1, 0, 1, 0, 0, 1};
    EXPECT_FALSE(eval_alpha(senses));
    EXPECT_TRUE(eval_beta(senses));
    EXPECT_TRUE(eval_gamma(senses));
    EXPECT_TRUE(eval_everywhere(senses));

    // Should evaluate to true (inside delta)
    senses = {0, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    EXPECT_TRUE(eval_delta(senses));
    EXPECT_TRUE(eval_everywhere(senses));
}
