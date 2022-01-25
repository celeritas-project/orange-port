//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstRandomProto.cc
 * \brief Tests for class RandomProto
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../RandomProto.hh"

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "celeritas_test.hh"
#include "orange/query/UnitMetadata.hh"
#include "orange/track/UnitTracker.hh"
#include "../CuboidShape.hh"
#include "../ShapeContainer.hh"
#include "ProtoTest.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class RandomProtoTest : public celeritas::ProtoTest
{
  protected:
    //! Box unit with a single particle type
    SPConstProto make_box_single() const;

    //! Box unit with multiple particle types
    SPConstProto make_box_multi() const;
};

//---------------------------------------------------------------------------//
/*!
 * Box unit with a single particle type
 */
auto RandomProtoTest::make_box_single() const -> SPConstProto
{
    RandomProto::Particle p1;
    p1.proto = this->make_sphere_unit(ORANGE_MD_FROM_SOURCE("sphere"),
                                      0.6203504908994);
    p1.volume_fraction = 0.1;
    p1.md              = ORANGE_MD_FROM_SOURCE("p1");

    // Make box big enough to fit a unit sphere with 30% packing
    ShapeContainer shapes;
    auto box = shapes.emplace<CuboidShape>(ORANGE_MD_FROM_SOURCE("cuboid"),
                                           Transform{},
                                           Real3{0, 0, 0},
                                           Real3{2, 5, 10});

    RandomProto::Params inp;
    inp.interior   = {{inside, box}};
    inp.fill_matid = 123; // arbitrary
    inp.particles  = {std::move(p1)};
    inp.md         = ORANGE_MD_FROM_SOURCE("rando_single");

    return std::make_shared<RandomProto>(std::move(inp));
}

//---------------------------------------------------------------------------//
/*!
 * Box unit with multiple particle types
 */
auto RandomProtoTest::make_box_multi() const -> SPConstProto
{
    // sphere w/v=1: 5% packed should produce 5
    RandomProto::Particle p1;
    p1.proto = this->make_sphere_unit(ORANGE_MD_FROM_SOURCE("sphere1"),
                                      0.6203504908994);
    p1.volume_fraction = 0.05;
    p1.md              = ORANGE_MD_FROM_SOURCE("p1");

    // sphere w/v=2 10% packed should produce 5
    RandomProto::Particle p2;
    p2.proto = this->make_sphere_unit(ORANGE_MD_FROM_SOURCE("sphere2"),
                                      0.781592641796772);
    p2.volume_fraction = 0.1;
    p2.md              = ORANGE_MD_FROM_SOURCE("p2");

    // Make box big enough to fit a unit sphere with 30% packing
    ShapeContainer shapes;
    auto box = shapes.emplace<CuboidShape>(ORANGE_MD_FROM_SOURCE("cuboid"),
                                           Transform{},
                                           Real3{0, 0, 0},
                                           Real3{2, 5, 10});

    RandomProto::Params inp;
    inp.interior   = {{inside, box}};
    inp.fill_matid = 1; // arbitrary
    inp.particles  = {std::move(p1), std::move(p2)};
    inp.md         = ORANGE_MD_FROM_SOURCE("rando_multi");

    return std::make_shared<RandomProto>(std::move(inp));
}

//---------------------------------------------------------------------------//

struct UnpackedBuildResult
{
    std::shared_ptr<const UnitTracker>  tracker;
    std::shared_ptr<const UnitMetadata> md;

    std::vector<std::string> fill;

    void print_expected() const;
};

void UnpackedBuildResult::print_expected() const
{
    using std::cout;

    cout << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
         << "static const char* const expected_fill[] = "
         << to_string(this->fill) << ";\n"
         << "EXPECT_VEC_EQ(expected_fill, built.fill);\n"
         << "/*** END CODE ***/\n";

    if (this->md && this->tracker)
    {
        cout << to_stream(*this->md, *this->tracker);
    }
    else
    {
        cout << "// WARNING: missing tracker and/or md\n";
    }
}

UnpackedBuildResult unpack_proto_build(Proto::BuildResult build)
{
    UnpackedBuildResult result;
    {
        // Save tracker
        std::shared_ptr<const Tracker> tracker(std::move(build.tracker));
        result.tracker = std::dynamic_pointer_cast<const UnitTracker>(tracker);
        EXPECT_TRUE(result.tracker);
    }
    {
        // Save metadata
        result.md = std::dynamic_pointer_cast<const UnitMetadata>(build.md);
        EXPECT_TRUE(result.md);
    }
    {
        // Fills
        CELER_ASSERT(build.md);
        result.fill.assign(build.md->num_volumes(), "<empty>");
        for (const auto& cell_matid : build.matids)
        {
            unsigned int cell = cell_matid.first.unchecked_get();
            EXPECT_LT(cell, result.fill.size());
            result.fill[cell] = std::to_string(cell_matid.second);
        }
        for (const auto& cell_daughter : build.daughters)
        {
            unsigned int cell = cell_daughter.first.unchecked_get();
            EXPECT_LT(cell, result.fill.size());
            CELER_ASSERT(cell_daughter.second.proto);
            result.fill[cell] = cell_daughter.second.proto->metadata().name();
        }
    }
    return result;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(RandomProtoTest, random_single_type_geom)
{
    SPConstProto unit = this->make_box_single();
    EXPECT_EQ("rando_single", unit->metadata().name());

    const auto& interior = unit->interior();
    ASSERT_EQ(1, interior.size());
    EXPECT_EQ(inside, interior.front().first);
    EXPECT_EQ("cuboid", interior.front().second->name());

    // Build the tracker/MD from the proto
    Proto::BuildArgs args;
    args.implicit_boundary         = true; // we're inside another unit
    auto                     built = unpack_proto_build(unit->build(args));
    static const char* const expected_fill[] = {"<empty>",
                                                "sphere",
                                                "sphere",
                                                "sphere",
                                                "sphere",
                                                "sphere",
                                                "sphere",
                                                "sphere",
                                                "sphere",
                                                "sphere",
                                                "sphere",
                                                "123"};
    EXPECT_VEC_EQ(expected_fill, built.fill);
}

//---------------------------------------------------------------------------//

TEST_F(RandomProtoTest, random_multi_type_geom)
{
    SPConstProto unit = this->make_box_multi();
    EXPECT_EQ("rando_multi", unit->metadata().name());

    const auto& interior = unit->interior();
    ASSERT_EQ(1, interior.size());
    EXPECT_EQ(inside, interior.front().first);
    EXPECT_EQ("cuboid", interior.front().second->name());

    // Build the tracker/MD from the proto
    Proto::BuildArgs args;
    args.implicit_boundary         = true; // we're inside another unit
    auto                     built = unpack_proto_build(unit->build(args));
    static const char* const expected_fill[] = {"<empty>",
                                                "sphere2",
                                                "sphere2",
                                                "sphere2",
                                                "sphere2",
                                                "sphere2",
                                                "sphere1",
                                                "sphere1",
                                                "sphere1",
                                                "sphere1",
                                                "sphere1",
                                                "1"};
    EXPECT_VEC_EQ(expected_fill, built.fill);
}
