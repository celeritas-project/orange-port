//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstUnionShape.cc
 * \brief Tests for class UnionShape
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../UnionShape.hh"

#include "celeritas_test.hh"
#include "base/Provenance.hh"
#include "orange/TransformUtils.hh"
#include "orange/query/ObjectMetadata.hh"
#include "orange/construct/PlacedShape.hh"
#include "orange/construct/SphereShape.hh"
#include "orange/surfaces/Sphere.hh"
#include "ShapeTest.hh"

using celeritas::UnionShape;

using celeritas::neg;
using celeritas::ObjectMetadata;
using celeritas::PlacedShape;
using celeritas::pos;
using celeritas::SphereShape;
using geometria::rotation_matrix;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class UnionShapeTest : public ::celeritas::ShapeTest
{
  public:
    using RegionVec    = UnionShape::RegionVec;
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

TEST_F(UnionShapeTest, two_daughters)
{
    auto a_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("a"), Transform{Real3(1, 1, 0)}, 3.0);
    auto b_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("b"), Transform{Real3(-1, -1, 0)}, 3.0);

    UnionShape shape(RegionVec{{neg, a_sphere}, {neg, b_sphere}});

    const RegionVec& interior = shape.interior();
    ASSERT_EQ(2, interior.size());
    EXPECT_EQ(neg, interior[0].first);
    EXPECT_EQ(a_sphere, interior[0].second);
    EXPECT_EQ(neg, interior[1].first);
    EXPECT_EQ(b_sphere, interior[1].second);

    EXPECT_EQ("union", shape.type());
    EXPECT_SOFT_EQ(0.0, shape.volume());
    EXPECT_FALSE(shape.is_convex());
    EXPECT_SOFT_EQ(0.0, shape.inradius());

    this->build(shape);

    static const logic_int expected_logic[]
        = {0l, LOGIC_NOT, 1l, LOGIC_NOT, LOGIC_OR};
    EXPECT_VEC_EQ(expected_logic, this->cell.logic());

    EXPECT_VEC_SOFT_EQ(Real3({-4, -4, -3}), this->bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3({4, 4, 3}), this->bbox.upper());

    static const char* expected_surface_names[] = {"s", "s"};
    EXPECT_VEC_EQ(expected_surface_names, this->surface_names);

    ASSERT_EQ(2, this->surfaces.size());
    EXPECT_VEC_SOFT_EQ(VecDbl({1, 1, 0, 9}),
                       this->get_surface_data(SurfaceType::s, 0));
    EXPECT_VEC_SOFT_EQ(VecDbl({-1, -1, 0, 9}),
                       this->get_surface_data(SurfaceType::s, 1));
}

TEST_F(UnionShapeTest, translated)
{
    auto a_sphere = this->make_shape<SphereShape>(
        ORANGE_MD_FROM_SOURCE("sphere"), Transform{}, 3.0);

    UnionShape shape(RegionVec{{neg, a_sphere}});

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
