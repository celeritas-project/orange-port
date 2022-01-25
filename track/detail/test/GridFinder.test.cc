//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/test/tstGridFinder.cc
 * \brief Tests for class GridFinder
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../GridFinder.hh"

#include <vector>
#include "Nemesis/containers/Span.hh"
#include "celeritas_test.hh"

using make_span;
using celeritas::detail::GridFinder;

GridFinder<int>::result_type make_result(unsigned int cell, int edge)
{
    CELER_EXPECT(edge >= -1 && edge <= 1);
    return {cell, edge};
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(GridFinderTest, all)
{
    std::vector<int> values = {1, 2, 4, 8, 16};
    GridFinder<int>  find{make_span(values)};

    EXPECT_EQ(make_result(0, -1), find(-123));
    EXPECT_EQ(make_result(0, -1), find(1));
    EXPECT_EQ(make_result(1, -1), find(2));
    EXPECT_EQ(make_result(1, 0), find(3));
    EXPECT_EQ(make_result(2, -1), find(4));
    EXPECT_EQ(make_result(3, 0), find(15));
    EXPECT_EQ(make_result(3, 1), find(16));
    EXPECT_EQ(make_result(3, 1), find(17));
}
