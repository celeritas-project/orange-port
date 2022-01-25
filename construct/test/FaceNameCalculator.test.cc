//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstFaceNameCalculator.cc
 * \brief Tests for class FaceNameCalculator
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../FaceNameCalculator.hh"

#include "celeritas_test.hh"
#include "../SlabShape.hh"
#include "../CuboidShape.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(FaceNameCalculatorTest, all)
{
    FaceNameCalculator calc_names;
    {
        auto names = calc_names(SlabShape(Axis::z, -3, 3));

        const std::string expected_names[] = {"mz", "pz"};
        EXPECT_VEC_EQ(expected_names, names);
    }

    {
        auto names = calc_names(CuboidShape({0, 0, 0}, {1, 2, 3}));

        const std::string expected_names[]
            = {"mx", "px", "my", "py", "mz", "pz"};
        EXPECT_VEC_EQ(expected_names, names);
    }
}
