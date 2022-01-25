//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstGeometryBuilder.cc
 * \brief Tests for class GeometryBuilder
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../GeometryBuilder.hh"

#include "celeritas_test.hh"
#include "orange/construct/SphereShape.hh"
#include "orange/query/ObjectMetadata.hh"
#include "../ShapeContainer.hh"
#include "../CuboidShape.hh"
#include "../SphereShape.hh"
#include "../UnitProto.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class GeometryBuilderTest : public ::Test
{
  protected:
    //// TYPE ALIASES ////
    using SPUnitProto = std::shared_ptr<UnitProto>;

  protected:
    SPUnitProto build_spheres_proto(real_type      radius,
                                    SPUnitProto    left,
                                    SPUnitProto    right,
                                    ObjectMetadata md) const
    {
        ShapeContainer shapes;

        // Build two spheres with scaled radius at the given position
        auto outer = shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("outer"), Transform{}, radius);

        Proto::RegionVec interior_def = {{inside, outer}};

        UnitProto::Params params;
        if (left)
        {
            UnitProto::Hole hole;
            hole.zorder    = ZOrder::media;
            hole.proto     = left;
            hole.transform = Real3{-radius / 2, 0, 0};
            hole.md        = ORANGE_MD_FROM_SOURCE("left hole");
            interior_def.push_back({pos, UnitProto::make_hole_shape(hole)});
            params.holes.push_back(std::move(hole));
        }
        if (right)
        {
            UnitProto::Hole hole;
            hole.zorder    = ZOrder::media;
            hole.proto     = right;
            hole.transform = Real3{radius / 2, 0, 0};
            hole.md        = ORANGE_MD_FROM_SOURCE("right hole");
            interior_def.push_back({pos, UnitProto::make_hole_shape(hole)});
            params.holes.push_back(std::move(hole));
        }
        {
            UnitProto::Media media;
            media.interior = std::move(interior_def);
            media.matid    = 1;
            media.md       = ORANGE_MD_FROM_SOURCE("sphere medium");
            params.media.push_back(std::move(media));
        }
        {
            UnitProto::Boundary boundary;
            boundary.interior          = {{inside, outer}};
            boundary.implicit_boundary = false; // connected to media
            boundary.md                = ORANGE_MD_FROM_SOURCE("boundary");
            params.boundary            = std::move(boundary);
        }
        params.md = std::move(md);
        return std::make_shared<UnitProto>(std::move(params));
    }

  protected:
    //// DATA ////
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(GeometryBuilderTest, simplest)
{
    auto global = build_spheres_proto(
        1.0, nullptr, nullptr, ORANGE_MD_FROM_SOURCE("global"));

    GeometryBuilder build_geo;
    auto            built = build_geo(global, {});

    EXPECT_EQ(1, built.geometry_params.universes.size());
    EXPECT_EQ(1, built.geometry_params.md.size());
    EXPECT_EQ(2, built.geometry_params.matids.size());
    EXPECT_EQ(0, built.geometry_params.boundaries.size());
    ASSERT_EQ(1, built.proto_to_univ.size());
    ASSERT_EQ(1, built.proto_to_univ.count(global));
    EXPECT_EQ(UniverseId{0}, built.proto_to_univ[global]);
}

TEST_F(GeometryBuilderTest, boundaries)
{
    GeometryBuilder::SPConstProto         global;
    GeometryBuilder::MapShapeFaceBoundary boundaries;

    {
        UnitProto::Params params;

        auto cuboid = std::make_shared<PlacedShape>(PlacedShape::Params{
            std::make_shared<CuboidShape>(Real3{0, 0, 0}, Real3{1, 2, 3}),
            {},
            ORANGE_MD_FROM_SOURCE("cub")});
        boundaries[{cuboid, "mx"}] = geometria::REFLECT;
        boundaries[{cuboid, "my"}] = geometria::REFLECT;

        {
            UnitProto::Media media;
            media.interior = {{inside, cuboid}};
            media.matid    = 1;
            media.md       = ORANGE_MD_FROM_SOURCE("cuboid medium");
            params.media.push_back(std::move(media));
        }
        {
            UnitProto::Boundary boundary;
            boundary.interior          = {{inside, cuboid}};
            boundary.implicit_boundary = true;
            boundary.md                = ORANGE_MD_FROM_SOURCE("boundary");
            params.boundary            = std::move(boundary);
        }
        params.md = ORANGE_MD_FROM_SOURCE("global");
        global    = std::make_shared<UnitProto>(std::move(params));
    }

    GeometryBuilder build_geo;
    auto            built = build_geo(global, boundaries);

    EXPECT_EQ(1, built.geometry_params.md.size());
    EXPECT_EQ(2, built.geometry_params.matids.size());

    const auto& bounds = built.geometry_params.boundaries;
    EXPECT_EQ(2, bounds.size());
    auto iter = built.geometry_params.boundaries.find(SurfaceId{0});
    ASSERT_NE(iter, built.geometry_params.boundaries.end());
    EXPECT_EQ(geometria::REFLECT, iter->second);
    iter = built.geometry_params.boundaries.find(SurfaceId{2});
    ASSERT_NE(iter, built.geometry_params.boundaries.end());
    EXPECT_EQ(geometria::REFLECT, iter->second);
}

TEST_F(GeometryBuilderTest, layered)
{
    auto micro = build_spheres_proto(
        1.0, nullptr, nullptr, ORANGE_MD_FROM_SOURCE("micro"));
    auto milli_dual = build_spheres_proto(
        2.0, micro, micro, ORANGE_MD_FROM_SOURCE("milli_dual"));
    auto milli_right = build_spheres_proto(
        2.0, nullptr, micro, ORANGE_MD_FROM_SOURCE("milli_right"));
    auto centi_combo = build_spheres_proto(
        4.0, milli_dual, milli_right, ORANGE_MD_FROM_SOURCE("centi_combo"));
    auto centi_left = build_spheres_proto(
        4.0, micro, nullptr, ORANGE_MD_FROM_SOURCE("centi_left"));
    auto global = build_spheres_proto(
        4.0, centi_left, centi_combo, ORANGE_MD_FROM_SOURCE("global"));

    GeometryBuilder build_geo;
    auto            built = build_geo(global, {});

    EXPECT_EQ(6, built.proto_to_univ.size());
    auto find = [&built](const SPUnitProto& p) -> UniverseId {
        auto iter = built.proto_to_univ.find(p);
        return iter == built.proto_to_univ.end() ? UniverseId{} : iter->second;
    };

    // IDs should be breadth-first based on the order of protos that we gave
    EXPECT_EQ(UniverseId{0}, find(global));
    EXPECT_EQ(UniverseId{1}, find(centi_left));
    EXPECT_EQ(UniverseId{2}, find(centi_combo));
    EXPECT_EQ(UniverseId{3}, find(micro));
    EXPECT_EQ(UniverseId{4}, find(milli_dual));
    EXPECT_EQ(UniverseId{5}, find(milli_right));
}
