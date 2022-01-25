//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/detail/test/tstSolveQuadratic.cc
 * \brief Utility function tests
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../SolveQuadratic.hh"

#include <cmath>
#include "celeritas_test.hh"

using celeritas::detail::Quadratic_Solver;
using celeritas::no_intersection();

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(SolveQuadraticTest, no_roots)
{
    // x^2 + 2*x + 1     = 0 -> one real root (x = -1)
    // x^2 + 2*x + 1.001 = 0 -> two complex roots
    {
        double b_2 = 1.0; // b/2
        double c   = 1.001;

        Quadratic_Solver solve_quadratic(1, b_2);
        auto             x = solve_quadratic(c);

        // Verify that x was not changed
        EXPECT_SOFT_EQ(no_intersection(), x[0]);
        EXPECT_SOFT_EQ(no_intersection(), x[1]);
    }

    // x^2 - 4*x + 4     = 0 -> one real root (x = 2)
    // x^2 - 4*x + 4.001 = 0 -> two complex roots
    {
        double b_2 = -2.0;
        double c   = 4.001;

        Quadratic_Solver solve_quadratic(1, b_2);
        auto             x = solve_quadratic(c);

        EXPECT_SOFT_EQ(no_intersection(), x[0]);
        EXPECT_SOFT_EQ(no_intersection(), x[1]);
    }
}

TEST(SolveQuadraticTest, one_root)
{
    // x^2 - 2*x + 1 = 0 -> x = 1
    double b_2 = -1.0;
    double c   = 1.0;

    // For solve_quadratic() to detect the single root case, it must calculate
    // that b/2 * b/2 - c == 0. It is assumed here that the floating point
    // operations -1.*-1. - 1. will reliably yield 0.

    Quadratic_Solver solve_quadratic(1, b_2);
    auto             x = solve_quadratic(c);

    EXPECT_SOFT_EQ(1.0, x[0]);
    EXPECT_SOFT_EQ(no_intersection(), x[1]);
}

TEST(SolveQuadraticTest, two_roots)
{
    // x^2 - 3*x + 2 = 0 -> x = 1, 2
    {
        double b_2 = -1.5;
        double c   = 2.0;

        auto x = Quadratic_Solver(1, b_2)(c);

        EXPECT_SOFT_EQ(1.0, x[0]);
        EXPECT_SOFT_EQ(2.0, x[1]);
    }

    // x^2 + x - 20 = 0 -> x = -5, 4
    {
        double b_2 = 0.5;
        double c   = -20.0;

        auto x = Quadratic_Solver(1, b_2)(c);

        EXPECT_SOFT_EQ(4.0, x[0]);
        EXPECT_SOFT_EQ(no_intersection(), x[1]);
    }

    // x^2 + 3*x + 2 = 0 -> x = -1, -2
    {
        double b_2 = 3.0 / 2;
        double c   = 2.0;

        auto x = Quadratic_Solver(1, b_2)(c);

        EXPECT_SOFT_EQ(no_intersection(), x[0]);
        EXPECT_SOFT_EQ(no_intersection(), x[1]);
    }

    // x^2 - 99999.999*x - 100 = 0 -> x = -0.001, 100000
    {
        double b_2 = -99999.999 / 2.0;
        double c   = -100.0;

        auto x = Quadratic_Solver(1, b_2)(c);

        EXPECT_SOFT_EQ(100000.0, x[0]);
        EXPECT_SOFT_EQ(no_intersection(), x[1]);
    }
}

//---------------------------------------------------------------------------//

TEST(SolveDegenerateQuadratic1Test, one_root)
{
    // x^2 + 2*x + 0 = 0, x -> -2
    {
        double b_2 = 2.0 / 2;

        auto x = Quadratic_Solver(1, b_2)();

        EXPECT_SOFT_EQ(no_intersection(), x[0]);
        EXPECT_SOFT_EQ(no_intersection(), x[1]);
    }

    // x^2 - 3.45*x + 0 = 0, x -> 3.45
    {
        double b_2 = -3.45 / 2;

        auto x = Quadratic_Solver(1, b_2)();

        EXPECT_SOFT_EQ(3.45, x[0]);
        EXPECT_SOFT_EQ(no_intersection(), x[1]);
    }
}

//---------------------------------------------------------------------------//

TEST(SolveDegenerateQuadratic2Test, no_roots)
{
    // For a*x^2 + b*x + c = 0, as both a and b -> 0, |x| -> infinity. In this
    // case, solve_degenerate_quadratic will return no positive roots.
    //
    // -1.0e-15*x^2 + 10000 = 0 -> x ~= +/- 3e9
    double a   = -1e-15;
    double b_2 = 0.;    // (b/a)/2
    double c   = 10000; // c/a

    auto x = Quadratic_Solver(a, b_2)(c);

    EXPECT_SOFT_EQ(no_intersection(), x[0]);
    EXPECT_SOFT_EQ(no_intersection(), x[1]);
}

TEST(SolveDegenerateQuadratic2Test, one_root)
{
    // 1.0e-15*x^2 + 2*x - 2000 = 0 -> x = 1e3
    double a   = 1e-15;
    double b_2 = 2.0 / 2;
    double c   = -2000 / a;

    auto x = Quadratic_Solver(a, b_2)(c);

    EXPECT_SOFT_EQ(1.0e3, x[0]);
    EXPECT_SOFT_EQ(no_intersection(), x[1]);
}
