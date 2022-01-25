//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/test/tstTransform.cc
 * \brief Tests for class Transform
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Transform.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "base/FixedViewArray.hh"
#include "base/VectorFunctions.hh"
#include "orange/TransformUtils.hh"
#include "orange/Transform.hh"

using make_fixed_view;
using celeritas::detail::transform_to_daughter;
using constants::two_pi;
using geometria::Real3;
using geometria::Transform;

//---------------------------------------------------------------------------//
// TRANSFORM_TO_DAUGHTER
//---------------------------------------------------------------------------//

TEST(TransformTest, to_daughter)
{
    Transform t(Real3{1.1, -half, 3.2});

    Real3 parent_pos{1, 20, 300};
    Real3 parent_dir{1, 2, 3};
    normalize_direction(&parent_dir);

    Real3 local_pos;
    Real3 local_dir;

    // Apply it
    transform_to_daughter(make_fixed_view(parent_pos),
                          make_fixed_view(parent_dir),
                          t,
                          make_fixed_view(local_pos),
                          make_fixed_view(local_dir));

    // Check
    const real_type expected_pos[] = {-0.1, 20.5, 296.8};

    EXPECT_VEC_SOFT_EQ(expected_pos, local_pos);
    EXPECT_VEC_SOFT_EQ(parent_dir, local_dir);
}

//---------------------------------------------------------------------------//

TEST(TransformTest, to_daughter_rotate)
{
    // Build a transform
    Transform t(
        geometria::rotation_matrix(
            Axis::z, 1. / 4, Axis::y, std::acos(0.2) / two_pi, Axis::x, 1. / 3),
        Real3(1.1, -half, 3.2));
    ASSERT_TRUE(t.has_rotation());

    Real3 parent_pos{1.071918358845, -0.6511911208929, -0.813272064836};
    Real3 parent_dir{0.6786799080846, -0.7260610035332, -0.1105848159138};
    Real3 local_pos;
    Real3 local_dir;

    // Apply it
    transform_to_daughter(make_fixed_view(parent_pos),
                          make_fixed_view(parent_dir),
                          t,
                          make_fixed_view(local_pos),
                          make_fixed_view(local_dir));

    // Check
    const real_type expected_pos[] = {-3.4, 2.1, 0.4};
    const real_type expected_dir[]
        = {0.2672612419124, 0.5345224838248, 0.8017837257373};

    EXPECT_VEC_NEAR(expected_pos, local_pos, 1e-10);
    EXPECT_VEC_NEAR(expected_dir, local_dir, 1e-10);
}
