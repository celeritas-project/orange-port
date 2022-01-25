//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstIntersectionShape.cc
 * \brief Tests for class IntersectionShape
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../IntersectionShape.hh"

#include "celeritas_test.hh"
#include "base/Provenance.hh"
#include "orange/TransformUtils.hh"
#include "orange/query/ObjectMetadata.hh"
#include "orange/construct/PlacedShape.hh"
#include "orange/construct/SphereShape.hh"
#include "orange/surfaces/Sphere.hh"
#include "ShapeTest.hh"

using celeritas::IntersectionShape;

using celeritas::neg;
using celeritas::ObjectMetadata;
using celeritas::PlacedShape;
using celeritas::pos;
using celeritas::SphereShape;
using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class IntersectionShapeTest : public ::celeritas::ShapeTest
{
  public:
    using RegionVec    = IntersectionShape::RegionVec;
    using SPConstShape = std::shared_ptr<const PlacedShape>;

    template<class S, typename... Args>
    SPConstShape
    make_shape(ObjectMetadata md, Transform transform, Args... args)
    {
        // Create shape
        auto shape = std::make_shared<S>(std::forward<Args>(args)...);

        // Create placed shape
        PlacedShape::Params params
            = {std::move(shape), std::move(transform), std::move(md)};

        return std::make_shared<PlacedShape>(std::move(params));
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(IntersectionShapeTest, two_daughters)
{
    auto a_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("a"), Transform{Real3(1, 1, 0)}, 3.0);
    auto b_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("b"), Transform{Real3(-1, -1, 0)}, 3.0);

    IntersectionShape shape(RegionVec{{neg, a_sphere}, {neg, b_sphere}});

    const RegionVec& interior = shape.interior();
    ASSERT_EQ(2, interior.size());
    EXPECT_EQ(neg, interior[0].first);
    EXPECT_EQ(a_sphere, interior[0].second);
    EXPECT_EQ(neg, interior[1].first);
    EXPECT_EQ(b_sphere, interior[1].second);

    EXPECT_EQ("intersection", shape.type());
    EXPECT_SOFT_EQ(0.0, shape.volume());
    EXPECT_TRUE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_NOT, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -2, -3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2, 2, 3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"s", "s"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(2, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 1, 0, 9}),
                       this->get_surface_data(SurfaceType::s, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1, -1, 0, 9}),
                       this->get_surface_data(SurfaceType::s, 1));
}

TEST_F(IntersectionShapeTest, translated)
{
    auto a_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("sphere"), Transform{}, 3.0);

    IntersectionShape shape(RegionVec{{neg, a_sphere}});

    this->build(shape, Transform{{1, 1, 0}});

    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -2, -3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4, 4, 3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"s"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 1, 0, 9}),
                       this->get_surface_data(SurfaceType::s, 0));
}

TEST_F(IntersectionShapeTest, neg_pos)
{
    auto a_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("a"), Transform{Real3(1, 1, 0)}, 3.0);
    auto b_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("b"), Transform{Real3(-1, -1, 0)}, 3.0);

    IntersectionShape shape(RegionVec{{neg, a_sphere}, {pos, b_sphere}});

    this->build(shape);

    static const logic_int expected_logic[] = {0l, LOGIC_NOT, 1l, LOGIC_AND};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -2, -3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4, 4, 3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"s", "s"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(2, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 1, 0, 9}),
                       this->get_surface_data(SurfaceType::s, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1, -1, 0, 9}),
                       this->get_surface_data(SurfaceType::s, 1));
}

TEST_F(IntersectionShapeTest, we_need_to_go_deeper)
{
    auto a_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("sphere"), Transform{}, 3.0);

    auto placed_meta
        = this->make_shape<IntersectionShape>(ORANGE_MD_FROM_SOURCE("meta"),
                                              Transform{Real3(1, 1, 0)},
                                              RegionVec{{neg, a_sphere}});

    // Create a meta-meta-shape
    auto inception = IntersectionShape(RegionVec{{neg, placed_meta}});

    this->build(inception,
                Transform{rotation_matrix(Axis::x, .25), {0, 0, 10}});

    static const logic_int expected_logic[] = {0l, LOGIC_NOT};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-2, -3, 8}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4, 3, 14}), this->bbox.upper());

    static const char* expected_surface_names[] = {"s"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(1, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 0, 11, 9}),
                       this->get_surface_data(SurfaceType::s, 0));
}
