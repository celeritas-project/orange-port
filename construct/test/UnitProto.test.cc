//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstUnitProto.cc
 * \brief Tests for class UnitProto
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../UnitProto.hh"

#include "celeritas_test.hh"
#include "base/StringFunctions.hh"
#include "orange/construct/UnitRegion.hh"
#include "orange/surfaces/CenteredSphere.hh"
#include "orange/surfaces/Sphere.hh"
#include "orange/query/UnitMetadata.hh"
#include "orange/track/UnitTracker.hh"
#include "orange/track/SimpleUnitTracker.hh"
#include "orange/track/MaskedUnitTracker.hh"
#include "../CuboidShape.hh"
#include "../IntersectionShape.hh"
#include "../ShapeContainer.hh"
#include "../SphereShape.hh"
#include "ProtoTest.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class UnitProtoTest : public celeritas::ProtoTest
{
  protected:
    // Build a proto with two spheres
    SPConstProto make_apple_unit() const;

    // Build a cuboid that simulates an array from (0, 0, 0) to (4, 5, 6)
    SPConstProto make_fake_array() const;
};

//---------------------------------------------------------------------------//
/*!
 * Build a unit that looks like an apple.
 *
 * The main shape is a sphere of radius three, and there's a small sphere that
 * overlaps on the boundary.
 */
auto UnitProtoTest::make_apple_unit() const -> SPConstProto
{
    ShapeContainer shapes;

    auto apple = shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("large sphere"), Transform{}, 3.0);
    auto bite = shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("bite sphere"), Transform{{3, 3, 0}}, 1.0);

    UnitProto::Params params;
    {
        UnitProto::Media media;
        media.interior = {{inside, apple}, {outside, bite}};
        media.matid    = 2;
        media.md       = ORANGE_MD_FROM_SOURCE("apple medium");
        params.media.push_back(std::move(media));
    }
    {
        UnitProto::Boundary boundary;
        boundary.interior          = {{inside, apple}, {outside, bite}};
        boundary.implicit_boundary = false; // connected to media
        boundary.md     = ORANGE_MD_FROM_SOURCE("sphere unit boundary");
        params.boundary = boundary;
    }
    params.md = ORANGE_MD_FROM_SOURCE("apple unit");
    return std::make_shared<UnitProto>(std::move(params));
}

//---------------------------------------------------------------------------//
/*!
 * Build a cuboid that simulates an array from (0, 0, 0) to (4, 5, 6).
 */
auto UnitProtoTest::make_fake_array() const -> SPConstProto
{
    ShapeContainer shapes;

    auto box = shapes.emplace<CuboidShape>(ORANGE_MD_FROM_SOURCE("cuboid"),
                                           Transform{},
                                           Real3{0, 0, 0},
                                           Real3{4, 5, 6});
    UnitProto::Params params;
    {
        UnitProto::Media media;
        media.interior = {{inside, box}};
        media.volume   = box->shape()->volume();
        media.md       = ORANGE_MD_FROM_SOURCE("box interior");
        media.matid    = 1234;
        params.media.push_back(std::move(media));
    }
    {
        UnitProto::Boundary boundary;
        boundary.interior          = {{inside, box}};
        boundary.implicit_boundary = false; // connected to media
        boundary.md     = ORANGE_MD_FROM_SOURCE("array unit boundary");
        params.boundary = boundary;
    }
    params.md = ORANGE_MD_FROM_SOURCE("array unit");
    return std::make_shared<UnitProto>(std::move(params));
}

//---------------------------------------------------------------------------//

struct UnpackedBuildResult
{
    std::shared_ptr<const UnitTracker>  tracker;
    std::shared_ptr<const UnitMetadata> md;

    std::vector<std::string> fill;
    std::vector<std::string> cells;
    std::vector<std::string> surfaces;

    void print_expected() const;
};

void UnpackedBuildResult::print_expected() const
{
    using std::cout;

    cout << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
         << "static const char* const expected_fill[] = "
         << to_string(this->fill) << ";\n"
         << "EXPECT_VEC_EQ(expected_fill, built.fill);\n"
         << "static const char* const expected_cells[] = "
         << to_string(this->cells) << ";\n"
         << "EXPECT_VEC_EQ(expected_cells, built.cells);\n"
         << "static const char* const expected_surfaces[] = "
         << to_string(this->surfaces) << ";\n"
         << "EXPECT_VEC_EQ(expected_surfaces, built.surfaces);\n"
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
    for (auto c : range(result.md->num_volumes()))
    {
        result.cells.push_back(result.md->id_to_label(VolumeId{c}));
    }
    for (auto s : range(result.md->num_surfaces()))
    {
        result.surfaces.push_back(result.md->id_to_label(SurfaceId{s}));
    }
    return result;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(UnitProtoTest, accessors)
{
    SPConstProto unit = this->make_sphere_unit(ORANGE_MD_FROM_SOURCE("su"));
    EXPECT_EQ("su", unit->metadata().name());

    const auto& interior = unit->interior();
    ASSERT_EQ(1, interior.size());
    EXPECT_EQ(inside, interior.front().first);
    EXPECT_EQ("sphere", interior.front().second->name());
}

//---------------------------------------------------------------------------//

TEST_F(UnitProtoTest, make_hole_shape)
{
    // Define daughter with *single* bounding surface in parent unit
    UnitProto::Hole hole;
    hole.proto     = this->make_sphere_unit(ORANGE_MD_FROM_SOURCE("su"));
    hole.transform = Real3{3, 0, 0}; // parent-to-daughter
    hole.md        = ORANGE_MD_FROM_SOURCE("sphere hole");

    auto hole_shape = UnitProto::make_hole_shape(hole);
    ASSERT_SP_TRUE(hole_shape);
    EXPECT_EQ("sphere hole", hole_shape->name());
    EXPECT_VEC_EQ(Real3(3, 0, 0), hole_shape->transform().translation());
    {
        auto raw_shape = std::dynamic_pointer_cast<const SphereShape>(
            hole_shape->shape());
        ASSERT_SP_TRUE(raw_shape);
        EXPECT_DOUBLE_EQ(3.0, raw_shape->radius());
    }

    // Define daughter with *two* bounding surfaces placed in parent unit
    hole.proto = this->make_apple_unit();
    hole.md    = ORANGE_MD_FROM_SOURCE("apple hole");

    hole_shape = UnitProto::make_hole_shape(hole);
    ASSERT_SP_TRUE(hole_shape);
    EXPECT_EQ("apple hole", hole_shape->name());
    EXPECT_EQ(Real3(3, 0, 0), hole_shape->transform().translation());
    {
        auto raw_shape = std::dynamic_pointer_cast<const IntersectionShape>(
            hole_shape->shape());
        ASSERT_SP_TRUE(raw_shape);
        const auto& interior = raw_shape->interior();
        ASSERT_EQ(2, interior.size());
        EXPECT_EQ(inside, interior[0].first);
        EXPECT_EQ("large sphere", interior[0].second->name());

        EXPECT_EQ(outside, interior[1].first);
        EXPECT_EQ("bite sphere", interior[1].second->name());
    }
}

//---------------------------------------------------------------------------//
// Build a well-connected global unit with a two-shape boundary
TEST_F(UnitProtoTest, build_global)
{
    auto proto = this->make_apple_unit();

    // Build the tracker/MD from the proto
    Proto::BuildArgs args;
    args.implicit_boundary = false; // always false for global unit
    auto built             = unpack_proto_build(proto->build(args));

    // built.print_expected();
    static const char* const expected_fill[] = {"<empty>", "2"};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    static const char* const expected_cells[]
        = {"sphere unit boundary", "apple medium"};
    EXPECT_VEC_EQ(expected_cells, built.cells);
    static const char* const expected_surfaces[]
        = {"large sphere.s", "bite sphere.s"};

    // Tracker should be "simple" since it's well connected and the boundary is
    // explicit
    ASSERT_SP_TRUE(built.tracker);
    const auto* tracker
        = dynamic_cast<const SimpleUnitTracker*>(built.tracker.get());
    ASSERT_SP_TRUE(tracker);
}

//---------------------------------------------------------------------------//
/*!
 * Build a global unit with a single daughter.
 */
TEST_F(UnitProtoTest, build_global_hole_parent)
{
    SPConstProto parent;
    {
        ShapeContainer shapes;

        auto sph = shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("parent sphere"), Transform{}, 5.0);

        UnitProto::Params params;
        {
            UnitProto::Media media;
            media.interior = {{inside, sph}};
            media.matid    = 4;
            media.volume   = sph->shape()->volume();
            media.md       = ORANGE_MD_FROM_SOURCE("parent medium");
            params.media.push_back(std::move(media));
        }
        {
            UnitProto::Hole hole;
            hole.proto
                = this->make_sphere_unit(ORANGE_MD_FROM_SOURCE("shole"));
            hole.transform = Real3{2, 0, 0};
            hole.md        = ORANGE_MD_FROM_SOURCE("translated hole");
            params.holes.push_back(std::move(hole));
        }
        {
            UnitProto::Boundary boundary;
            boundary.interior          = {{inside, sph}};
            boundary.implicit_boundary = true;
            boundary.md     = ORANGE_MD_FROM_SOURCE("parent unit boundary");
            params.boundary = boundary;
        }
        params.md = ORANGE_MD_FROM_SOURCE("parent unit");
        parent    = std::make_shared<UnitProto>(std::move(params));
    }
    CELER_ASSERT(parent);

    // Build the tracker/MD from the proto
    Proto::BuildArgs args;
    args.implicit_boundary = false; // always false for global unit
    auto built             = unpack_proto_build(parent->build(args));

    // built.print_expected();
    static const char* const expected_fill[] = {"<empty>", "shole", "4"};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    static const char* const expected_cells[]
        = {"parent unit boundary", "translated hole", "parent medium"};
    EXPECT_VEC_EQ(expected_cells, built.cells);
    static const char* const expected_surfaces[]
        = {"parent sphere.s", "translated hole.s"};
    EXPECT_VEC_EQ(expected_surfaces, built.surfaces);

    // Tracker should be "masked" since we're using holes that implicitly
    // truncate the media
    ASSERT_SP_TRUE(built.tracker);
    const auto* tracker
        = dynamic_cast<const MaskedUnitTracker*>(built.tracker.get());
    ASSERT_SP_TRUE(tracker);

    // Check regions
    {
        // Exterior region
        auto reg = tracker->get_region(VolumeId{0});
        EXPECT_EQ(CSGCell::from_string("0").logic(), reg.interior.logic());
        EXPECT_EQ(ZOrder::exterior, reg.zorder);
    }
    {
        // Hole region
        auto reg = tracker->get_region(VolumeId{1});
        EXPECT_EQ(CSGCell::from_string("1 ~").logic(), reg.interior.logic());
        EXPECT_EQ(ZOrder::hole, reg.zorder);
    }
    {
        // Fill
        auto reg = tracker->get_region(VolumeId{2});
        EXPECT_EQ(CSGCell::from_string("0 ~").logic(), reg.interior.logic());
        EXPECT_EQ(ZOrder::media, reg.zorder);
    }

    // Check surfaces
    const auto& surfaces = tracker->surfaces();
    {
        EXPECT_EQ(SurfaceType::so, surfaces.get_type(SurfaceId{0}));
        auto surf = surfaces.get<CenteredSphere>(SurfaceId{0});
        EXPECT_SOFT_EQ(5.0 * 5.0, surf.radius_sq());
    }
    {
        EXPECT_EQ(SurfaceType::s, surfaces.get_type(SurfaceId{1}));
        auto surf = surfaces.get<Sphere>(SurfaceId{1});
        EXPECT_SOFT_EQ(3.0 * 3.0, surf.radius_sq());
        EXPECT_VEC_SOFT_EQ(Real3(2, 0, 0), surf.origin());
    }
}

//---------------------------------------------------------------------------//
/*!
 * Build an interior unit with a single KENO-VI array daughter and a media.
 *
 * The "interior" shapes should be deleted from the media definitions since
 * they're *implicitly* imposed from the higher-order unit.
 */
TEST_F(UnitProtoTest, build_array_parent)
{
    SPConstProto parent;
    {
        ShapeContainer shapes;

        auto sph = shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("parent sphere"), Transform{}, 10.0);
        auto arrhole = shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("array placement"),
            Transform({2, 2.5, 3}),
            5.0);

        UnitProto::Params params;
        {
            UnitProto::Media media;
            media.interior = {{inside, sph}, {outside, arrhole}};
            media.matid    = 4;
            media.volume   = sph->shape()->volume();
            media.md       = ORANGE_MD_FROM_SOURCE("parent medium");
            params.media.push_back(std::move(media));
        }
        {
            UnitProto::Array array;
            array.proto     = this->make_fake_array();
            array.transform = Real3{0, 0, 0};
            array.md        = ORANGE_MD_FROM_SOURCE("array placement");
            array.interior  = {{inside, arrhole}};
            params.arrays.push_back(std::move(array));
        }
        {
            UnitProto::Boundary boundary;
            boundary.interior          = {{inside, sph}};
            boundary.implicit_boundary = true;
            boundary.md     = ORANGE_MD_FROM_SOURCE("parent unit boundary");
            params.boundary = boundary;
        }
        params.md = ORANGE_MD_FROM_SOURCE("parent unit");
        parent    = std::make_shared<UnitProto>(std::move(params));
    }
    CELER_ASSERT(parent);

    // Build the tracker/MD from the proto
    Proto::BuildArgs args;
    args.implicit_boundary = true; // we're inside another unit
    auto built             = unpack_proto_build(parent->build(args));

    ASSERT_SP_TRUE(built.md);
    EXPECT_EQ("parent unit", built.md->metadata().name());

    // cout << to_stream(built.md, *built.tracker);

    // Tracker should be "simple" since there is only a single level (all
    // media)
    ASSERT_SP_TRUE(built.tracker);
    const auto* tracker
        = dynamic_cast<const SimpleUnitTracker*>(built.tracker.get());
    ASSERT_SP_TRUE(tracker);

    static const char* const expected_fill[] = {"<empty>", "array unit", "4"};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    static const char* const expected_cells[]
        = {"parent unit boundary", "array placement", "parent medium"};
    EXPECT_VEC_EQ(expected_cells, built.cells);
    static const char* const expected_surfaces[]
        = {"array placement.s", "parent sphere.s"};
    EXPECT_VEC_EQ(expected_surfaces, built.surfaces);

    // Check regions
    {
        // Implicit exterior
        auto reg = tracker->get_region(VolumeId{0});
        EXPECT_EQ(CSGCell::from_string("* ~").logic(), reg.interior.logic());
        EXPECT_EQ(ZOrder::implicit_exterior, reg.zorder);
    }
    {
        // Array region
        auto reg = tracker->get_region(VolumeId{1});
        EXPECT_EQ(CSGCell::from_string("0 ~").logic(), reg.interior.logic());
        EXPECT_EQ(ZOrder::media, reg.zorder);
    }
    {
        // Media outside array
        auto reg = tracker->get_region(VolumeId{2});
        EXPECT_EQ(CSGCell::from_string("0").logic(), reg.interior.logic());
        EXPECT_EQ(ZOrder::media, reg.zorder);
    }

    // Check surfaces
    const auto& surfaces = tracker->surfaces();
    {
        EXPECT_EQ(SurfaceType::s, surfaces.get_type(SurfaceId{0}));
        auto surf = surfaces.get<Sphere>(SurfaceId{0});
        EXPECT_SOFT_EQ(5.0 * 5.0, surf.radius_sq());
        EXPECT_VEC_SOFT_EQ(Real3(2, 2.5, 3), surf.origin());
    }
}
