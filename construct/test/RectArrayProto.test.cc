//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstRectArrayProto.cc
 * \brief Tests for class RectArrayProto
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../RectArrayProto.hh"

#include <memory>

#include "celeritas_test.hh"
#include "base/StringFunctions.hh"
#include "orange/query/RectArrayMetadata.hh"
#include "orange/track/RectArrayTracker.hh"
#include "../CuboidShape.hh"
#include "../PlacedShape.hh"
#include "../SphereShape.hh"
#include "../UnitProto.hh"
#include "ProtoTest.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class RectArrayProtoTest : public celeritas::ProtoTest
{
  protected:
    //// TYPE ALIASES ////
    using SPConstArray = std::shared_ptr<const ArrayProto>;
    using ProtoVec3    = ArrayProto::ProtoVec3;

  public:
    // Make a simple unit with an untransformed cuboid
    SPConstProto
    make_array_cell(ObjectMetadata unit_md, Real3 width, Real3 lower) const
    {
        PlacedShape::Params ps_params;
        ps_params.transform = {};
        ps_params.md        = ORANGE_MD_FROM_SOURCE("cuboid");
        ps_params.shape = std::make_shared<CuboidShape>(lower, lower + width);
        return this->make_simple_unit(std::move(unit_md), ps_params);
    }

    // Make an array cell with lower point at the origin
    SPConstProto make_array_cell(ObjectMetadata unit_md, Real3 width) const
    {
        return this->make_array_cell(std::move(unit_md), width, {0, 0, 0});
    }
};

struct UnpackedBuildResult
{
    std::shared_ptr<const RectArrayTracker>  tracker;
    std::shared_ptr<const RectArrayMetadata> md;

    std::vector<real_type>   x, y, z;     //!< Grid edges
    std::vector<std::string> fill;        //!< Fill names
    std::vector<real_type>   translation; //!< Flattened cell translations

    void print_expected() const;
};

void UnpackedBuildResult::print_expected() const
{
    cout << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
         << "static const char* const expected_fill[] = "
         << to_string(this->fill) << ";\n"
         << "static const real_type expected_x[] = " << to_string(this->x)
         << ";\n"
         << "static const real_type expected_y[] = " << to_string(this->y)
         << ";\n"
         << "static const real_type expected_z[] = " << to_string(this->z)
         << ";\n"
         << "static const real_type expected_translation[] = "
         << to_string(this->translation) << ";\n"
         << "EXPECT_VEC_EQ(expected_fill, built.fill);\n"
         << "EXPECT_VEC_EQ(expected_x, built.x);\n"
         << "EXPECT_VEC_EQ(expected_y, built.y);\n"
         << "EXPECT_VEC_EQ(expected_z, built.z);\n"
         << "EXPECT_VEC_EQ(expected_translation, built.translation);\n"
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
        result.tracker
            = std::dynamic_pointer_cast<const RectArrayTracker>(tracker);
        EXPECT_TRUE(result.tracker);
    }
    {
        // Save metadata
        result.md
            = std::dynamic_pointer_cast<const RectArrayMetadata>(build.md);
        EXPECT_TRUE(result.md);
    }
    {
        // Fills
        CELER_ASSERT(build.md);
        result.fill.resize(build.md->num_volumes());
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

            const Transform& transform = cell_daughter.second.transform;
            EXPECT_FALSE(transform.has_rotation());
            const Real3& translation = transform.translation();
            result.translation.insert(result.translation.end(),
                                      translation.begin(),
                                      translation.end());
        }
    }

    if (result.tracker && result.md)
    {
        const auto&             grid = result.tracker->grid();
        std::vector<real_type>* result_grids[]
            = {&result.x, &result.y, &result.z};
        std::vector<real_type> temp_grid;
        auto                   bbox = result.md->bbox();

        for (auto ax : range(3))
        {
            // Outer points are bumped away but it doesn't matter if they're
            // infinite.
            temp_grid = grid.edges(ax);
            EXPECT_LT(temp_grid.front(), bbox.lower()[ax]);
            EXPECT_GT(temp_grid.back(), bbox.upper()[ax]);
            temp_grid.front() = bbox.lower()[ax];
            temp_grid.back()  = bbox.upper()[ax];
            *result_grids[ax] = std::move(temp_grid);
        }
    }

    return result;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(RectArrayProtoTest, single_unit)
{
    SPConstArray arr;
    {
        auto daughter
            = this->make_array_cell(ORANGE_MD_FROM_SOURCE("unit"), {1, 2, 3});
        RectArrayProto::Params params;
        params.units = ProtoVec3{{1, 1, 1}, ProtoVec3::Storage_t{daughter}};
        params.md    = ORANGE_MD_FROM_SOURCE("parent array");
        arr          = std::make_shared<RectArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    // Check proto metadata
    EXPECT_EQ("parent array", arr->metadata().name());

    // Check interior shape
    const auto& interior = arr->interior();
    ASSERT_EQ(1, interior.size());
    EXPECT_EQ(inside, interior.front().first);
    auto shape = std::dynamic_pointer_cast<const CuboidShape>(
        interior.front().second->shape());
    ASSERT_SP_TRUE(shape);
    EXPECT_VEC_SOFT_EQ(Real3({0, 0, 0}), shape->lower());
    EXPECT_VEC_SOFT_EQ(Real3({1, 2, 3}), shape->upper());

    // Check placement
    {
        // Lower-left of the array is at (10, 0, 0) in parent space
        auto transform = arr->calc_placement(Real3(10, 0, 0));
        EXPECT_FALSE(transform.has_rotation());
        EXPECT_VEC_SOFT_EQ(Real3(10, 0, 0), transform.translation());
    }
    {
        // Transform is from array cell center, not daughter
        auto transform = arr->calc_placement({0, 0, 0}, Real3{2, 0, 0});
        EXPECT_FALSE(transform.has_rotation());
        EXPECT_VEC_SOFT_EQ(Real3(2, 0, 0), transform.translation());
    }

    // Build
    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr->build(args));

    // Check fills etc.
    static const char* const expected_fill[]        = {"unit"};
    static const real_type   expected_x[]           = {0, 1};
    static const real_type   expected_y[]           = {0, 2};
    static const real_type   expected_z[]           = {0, 3};
    static const real_type   expected_translation[] = {0, 0, 0};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_EQ(expected_x, built.x);
    EXPECT_VEC_EQ(expected_y, built.y);
    EXPECT_VEC_EQ(expected_z, built.z);
    EXPECT_VEC_EQ(expected_translation, built.translation);

    // Check built metadata
    ASSERT_SP_TRUE(built.md);
    EXPECT_EQ("parent array", built.md->metadata().name());
}

//---------------------------------------------------------------------------//
//
//  Place the same unit (centered on origin) in 3x1x1 array
//
//         ------- ------- -------
//        |       |       |       |
//        |       |       |       |
//        |  11   |  11   |  11   |
//   ^ y  |       |       |       |
//   |     ------- ------- -------
//   |
//   '----> x
//
TEST_F(RectArrayProtoTest, placement_all_centered)
{
    // Cuboid centered on origin
    SPConstProto u11 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u11"), {20, 30, 40}, {-10, -15, -20});

    SPConstArray arr;
    {
        // 3x1x1 array filled with the same unit, u11
        RectArrayProto::Params params;
        params.units
            = ProtoVec3{{3, 1, 1}, ProtoVec3::Storage_t{u11, u11, u11}};
        params.md = ORANGE_MD_FROM_SOURCE("parent array");
        arr       = std::make_shared<RectArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);
    {
        auto calc_origin = [arr](RectArrayProto::DimVector ijk,
                                 SpanConstReal3            xyz) -> Real3 {
            Transform to_array = arr->calc_placement(ijk, xyz);
            return to_array.translation();
        };

        // Translation to 'cell centers' should all be zero
        Real3 origin{0, 0, 0};
        //
        EXPECT_VEC_EQ(origin, calc_origin({0, 0, 0}, {10.0, 15.0, 20.0}));
        EXPECT_VEC_EQ(origin, calc_origin({1, 0, 0}, {30.0, 15.0, 20.0}));
        EXPECT_VEC_EQ(origin, calc_origin({2, 0, 0}, {50.0, 15.0, 20.0}));
    }

    // Build
    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr->build(args));

    static const char* const expected_fill[] = {"u11", "u11", "u11"};
    static const real_type   expected_x[]    = {0, 20, 40, 60};
    static const real_type   expected_y[]    = {0, 30};
    static const real_type   expected_z[]    = {0, 40};
    static const real_type   expected_translation[]
        = {10, 15, 20, 30, 15, 20, 50, 15, 20};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_EQ(expected_x, built.x);
    EXPECT_VEC_EQ(expected_y, built.y);
    EXPECT_VEC_EQ(expected_z, built.z);
    EXPECT_VEC_EQ(expected_translation, built.translation);
}

//---------------------------------------------------------------------------//
//
//  Place the same unit (off-center) in 3x1x1 array
//
//         ------- ------- -------
//        |       |       |       |
//        |       |       |       |
//        |  11   |  11   |  11   |
//   ^ y  |       |       |       |
//   |     ------- ------- -------
//   |
//   '----> x
//
TEST_F(RectArrayProtoTest, placement_off_center)
{
    // Fill same cuboid (off-center)
    SPConstProto u11 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u11"), {20, 30, 40}, {100, 100, 100});

    SPConstArray arr;
    {
        // 3x1x1 array filled with the same unit, u11
        RectArrayProto::Params params;
        params.units
            = ProtoVec3{{3, 1, 1}, ProtoVec3::Storage_t{u11, u11, u11}};
        params.md = ORANGE_MD_FROM_SOURCE("parent array");
        arr       = std::make_shared<RectArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);
    {
        auto calc_origin = [arr](RectArrayProto::DimVector ijk,
                                 SpanConstReal3            xyz) -> Real3 {
            Transform to_array = arr->calc_placement(ijk, xyz);
            // Return translation
            return to_array.translation();
        };
        // Translation to 'cell centers' should all be zero
        Real3 origin{0, 0, 0};
        //
        EXPECT_VEC_EQ(origin, calc_origin({0, 0, 0}, {-100, -100, -100}));
        EXPECT_VEC_EQ(origin, calc_origin({1, 0, 0}, {-80, -100, -100}));
        EXPECT_VEC_EQ(origin, calc_origin({2, 0, 0}, {-60, -100, -100}));
    }

    // Build
    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr->build(args));

    static const char* const expected_fill[] = {"u11", "u11", "u11"};
    static const real_type   expected_x[]    = {0, 20, 40, 60};
    static const real_type   expected_y[]    = {0, 30};
    static const real_type   expected_z[]    = {0, 40};
    static const real_type   expected_translation[]
        = {-100, -100, -100, -80, -100, -100, -60, -100, -100};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_EQ(expected_x, built.x);
    EXPECT_VEC_EQ(expected_y, built.y);
    EXPECT_VEC_EQ(expected_z, built.z);
    EXPECT_VEC_EQ(expected_translation, built.translation);
}

//---------------------------------------------------------------------------//
//
// Place 6 units in 3x2x1 array (all centered on origin with different size
// in x and y)
//
//
//         --------- ----- ---------------
//        |         |     |               |
//        |         |     |               |
//        |  14     | 15  |    16         |
//        |         |     |               |
//        |         |     |               |
//        |         |     |               |
//        |         |     |               |
//         --------- ----- ---------------
//        |         |     |               |
//        |         |     |               |
//        |         |     |               |
//        |   11    | 12  |    13         |
//   ^ y  |         |     |               |
//   |    |         |     |               |
//   |     --------- ----- ---------------
//   ----> x

TEST_F(RectArrayProtoTest, placement_all_centered1)
{
    // Different size of unit cubes (in x and y) centered on origin
    SPConstProto u11 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u11"), {20, 30, 30}, {-10, -15, -15});
    SPConstProto u12 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u12"), {10, 30, 30}, {-5, -15, -15});
    SPConstProto u13 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u13"), {30, 30, 30}, {-15, -15, -15});
    //
    SPConstProto u14 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u14"), {20, 40, 30}, {-10, -20, -15});
    SPConstProto u15 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u15"), {10, 40, 30}, {-5, -20, -15});
    SPConstProto u16 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u16"), {30, 40, 30}, {-15, -20, -15});

    SPConstArray arr;
    {
        // 3x2x1 array
        RectArrayProto::Params params;
        params.units = ProtoVec3{
            {3, 2, 1}, ProtoVec3::Storage_t{u11, u14, u12, u15, u13, u16}};
        params.md = ORANGE_MD_FROM_SOURCE("parent array");
        arr       = std::make_shared<RectArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    {
        auto calc_origin = [arr](RectArrayProto::DimVector ijk,
                                 SpanConstReal3            xyz) -> Real3 {
            Transform to_array = arr->calc_placement(ijk, xyz);
            return to_array.translation();
        };

        // Translation to 'cell centers' should all be zero
        Real3 origin{0, 0, 0};
        //
        EXPECT_VEC_EQ(origin, calc_origin({0, 0, 0}, {10.0, 15.0, 15.0}));
        EXPECT_VEC_EQ(origin, calc_origin({1, 0, 0}, {25.0, 15.0, 15.0}));
        EXPECT_VEC_EQ(origin, calc_origin({2, 0, 0}, {45.0, 15.0, 15.0}));
        EXPECT_VEC_EQ(origin, calc_origin({0, 1, 0}, {10.0, 50.0, 15.0}));
        EXPECT_VEC_EQ(origin, calc_origin({1, 1, 0}, {25.0, 50.0, 15.0}));
        EXPECT_VEC_EQ(origin, calc_origin({2, 1, 0}, {45.0, 50.0, 15.0}));
    }

    // Build
    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr->build(args));

    static const char* const expected_fill[]
        = {"u11", "u12", "u13", "u14", "u15", "u16"};
    static const real_type expected_x[]           = {0, 20, 30, 60};
    static const real_type expected_y[]           = {0, 30, 70};
    static const real_type expected_z[]           = {0, 30};
    static const real_type expected_translation[] = {
        10, 15, 15, 25, 15, 15, 45, 15, 15, 10, 50, 15, 25, 50, 15, 45, 50, 15};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_EQ(expected_x, built.x);
    EXPECT_VEC_EQ(expected_y, built.y);
    EXPECT_VEC_EQ(expected_z, built.z);
    EXPECT_VEC_EQ(expected_translation, built.translation);
}

//---------------------------------------------------------------------------//
//
// Place 6 units in 3x2x1 array (all off-center with different size in x and y)
//
//
//         --------- ----- ---------------
//        |         |     |               |
//        |         |     |               |
//        |  14     | 15  |    16         |
//        |         |     |               |
//        |         |     |               |
//        |         |     |               |
//        |         |     |               |
//         --------- ----- ---------------
//        |         |     |               |
//        |         |     |               |
//        |         |     |               |
//        |   11    | 12  |    13         |
//   ^ y  |         |     |               |
//   |    |         |     |               |
//   |     --------- ----- ---------------
//   ----> x

TEST_F(RectArrayProtoTest, placement_off_center_2)
{
    // Different size of unit cubes in x and y, all off-center
    SPConstProto u11 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u11"), {20, 30, 30}, {-20, -115, 0});
    SPConstProto u12 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u12"), {10, 30, 30}, {-5, 0, -20});
    SPConstProto u13 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u13"), {30, 30, 30}, {-95, 120, 15});
    //
    SPConstProto u14 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u14"), {20, 40, 30}, {0, -120, -115});
    SPConstProto u15 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u15"), {10, 40, 30}, {-15, 120, 15});
    SPConstProto u16 = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("u16"), {30, 40, 30}, {15, 20, 100});

    SPConstArray arr;
    {
        // 3x2x1 array
        RectArrayProto::Params params;
        params.units = ProtoVec3{
            {3, 2, 1}, ProtoVec3::Storage_t{u11, u14, u12, u15, u13, u16}};
        params.md = ORANGE_MD_FROM_SOURCE("parent array");
        arr       = std::make_shared<RectArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    {
        auto calc_origin = [arr](RectArrayProto::DimVector ijk,
                                 SpanConstReal3            xyz) -> Real3 {
            Transform to_array = arr->calc_placement(ijk, xyz);
            return to_array.translation();
        };

        // Translation to 'cell centers' should all be zero
        Real3 origin{0, 0, 0};
        //
        EXPECT_VEC_EQ(origin, calc_origin({0, 0, 0}, {20, 115, 0}));
        EXPECT_VEC_EQ(origin, calc_origin({1, 0, 0}, {25, 0, 20}));
        EXPECT_VEC_EQ(origin, calc_origin({2, 0, 0}, {125, -120, -15}));
        EXPECT_VEC_EQ(origin, calc_origin({0, 1, 0}, {0, 150, 115}));
        EXPECT_VEC_EQ(origin, calc_origin({1, 1, 0}, {35, -90, -15}));
        EXPECT_VEC_EQ(origin, calc_origin({2, 1, 0}, {15, 10, -100}));
    }

    // Build
    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr->build(args));

    static const char* const expected_fill[]
        = {"u11", "u12", "u13", "u14", "u15", "u16"};
    static const real_type expected_x[]           = {0, 20, 30, 60};
    static const real_type expected_y[]           = {0, 30, 70};
    static const real_type expected_z[]           = {0, 30};
    static const real_type expected_translation[] = {20,
                                                     115,
                                                     0,
                                                     25,
                                                     0,
                                                     20,
                                                     125,
                                                     -120,
                                                     -15,
                                                     0,
                                                     150,
                                                     115,
                                                     35,
                                                     -90,
                                                     -15,
                                                     15,
                                                     10,
                                                     -100};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_EQ(expected_x, built.x);
    EXPECT_VEC_EQ(expected_y, built.y);
    EXPECT_VEC_EQ(expected_z, built.z);
    EXPECT_VEC_EQ(expected_translation, built.translation);
}

//---------------------------------------------------------------------------//

TEST_F(RectArrayProtoTest, oned_transforms)
{
    // Unit cube centered on origin
    SPConstProto centered_unit = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("centered"), {1, 1, 1}, {-.5, -.5, -.5});
    // Cube with length 2, right of origin
    SPConstProto right_unit = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("right"), {2, 1, 1}, {0, 0, 0});
    // Cube with length 3, left of origin
    SPConstProto left_unit = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("left"), {3, 1, 1}, {-3, 0, 0});
    // Transformed unit cube: lower is at x=10
    SPConstProto transformed_unit;
    {
        PlacedShape::Params ps_params;
        ps_params.transform = Transform(Real3{10, 0, 0});
        ps_params.md        = ORANGE_MD_FROM_SOURCE("cuboid");
        ps_params.shape
            = std::make_shared<CuboidShape>(Real3{0, 0, 0}, Real3{1, 1, 1});
        transformed_unit = this->make_simple_unit(
            ORANGE_MD_FROM_SOURCE("transformed"), std::move(ps_params));
    }

    SPConstArray arr;
    {
        // Unit with non-origin cuboid "lower" point
        RectArrayProto::Params params;
        params.units = ProtoVec3{
            {4, 1, 1},
            {centered_unit, right_unit, left_unit, transformed_unit}};
        params.md = ORANGE_MD_FROM_SOURCE("parent array");
        arr       = std::make_shared<RectArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    {
        auto calc_origin = [arr](unsigned int i, real_type x) -> real_type {
            Transform to_array = arr->calc_placement({i, 0, 0}, Real3{x, 0, 0});
            // Return X component of translation
            return to_array.translation()[0];
        };

        // Translation to 'cell centers' should all be zero
        EXPECT_SOFT_EQ(0.0, calc_origin(0, 0.5));
        EXPECT_SOFT_EQ(0.0, calc_origin(1, 1.0));
        EXPECT_SOFT_EQ(0.0, calc_origin(2, 6.0));
        EXPECT_SOFT_EQ(0.0, calc_origin(3, -4.0));

        EXPECT_SOFT_EQ(10.0, calc_origin(0, 10.5));
    }

    Proto::BuildArgs args;
    args.implicit_boundary         = true; // array is always truncated
    auto                     built = unpack_proto_build(arr->build(args));
    static const char* const expected_fill[]
        = {"centered", "right", "left", "transformed"};
    static const real_type expected_x[] = {0, 1, 3, 6, 7};
    static const real_type expected_y[] = {0, 1};
    static const real_type expected_z[] = {0, 1};
    static const real_type expected_translation[]
        = {0.5, 0.5, 0.5, 1, 0, 0, 6, 0, 0, -4, 0, 0};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_EQ(expected_x, built.x);
    EXPECT_VEC_EQ(expected_y, built.y);
    EXPECT_VEC_EQ(expected_z, built.z);
    EXPECT_VEC_EQ(expected_translation, built.translation);
}

//---------------------------------------------------------------------------//

TEST_F(RectArrayProtoTest, transforms)
{
    SPConstProto lower_unit = this->make_array_cell(
        ORANGE_MD_FROM_SOURCE("lower"), {1, 2, 3}, {1, 0, 0});

    SPConstProto transformed_unit;
    {
        // Boundary cuboid has a transform applied: its lower corner is
        // (0,0,-10)
        PlacedShape::Params ps_params;
        ps_params.transform = Transform(Real3{0, 0, -10});
        ps_params.md        = ORANGE_MD_FROM_SOURCE("cuboid");
        ps_params.shape
            = std::make_shared<CuboidShape>(Real3{0, 0, 0}, Real3{1, 4, 3});
        transformed_unit = this->make_simple_unit(
            ORANGE_MD_FROM_SOURCE("transformed"), std::move(ps_params));
    }

    SPConstArray arr;
    {
        // Unit with non-origin cuboid "lower" point
        RectArrayProto::Params params;
        params.units = ProtoVec3{{1, 2, 1}, {lower_unit, transformed_unit}};
        params.md    = ORANGE_MD_FROM_SOURCE("parent array");
        arr          = std::make_shared<RectArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr->build(args));

    static const char* const expected_fill[]        = {"lower", "transformed"};
    static const real_type   expected_x[]           = {0, 1};
    static const real_type   expected_y[]           = {0, 2, 6};
    static const real_type   expected_z[]           = {0, 3};
    static const real_type   expected_translation[] = {-1, 0, 0, 0, 2, 10};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_EQ(expected_x, built.x);
    EXPECT_VEC_EQ(expected_y, built.y);
    EXPECT_VEC_EQ(expected_z, built.z);
    EXPECT_VEC_EQ(expected_translation, built.translation);

    // Check transform in practice
    {
        // Place center of second unit at {3, 0, 0}.
        auto to_array = arr->calc_placement({0, 1, 0}, Real3{3, 0, 0});
        EXPECT_FALSE(to_array.has_rotation());
        EXPECT_VEC_SOFT_EQ(Real3(3.0, -2.0, -10.0), to_array.translation());

        // Build array
        auto built = arr->build(args);
        ASSERT_EQ(2, built.daughters.size());
        ASSERT_EQ(VolumeId{0}, built.daughters[0].first);
        ASSERT_EQ(VolumeId{1}, built.daughters[1].first);

        // Check that in practice, {3, 0, 0} in the global frame actually gets
        // transformed to the daughter's center
        Real3 pos = {3, 0, 0};
        to_array.parent_to_daughter(pos);
        built.daughters[1].second.transform.parent_to_daughter(pos);
        EXPECT_VEC_SOFT_EQ(Real3(0.0, 0.0, 0.0), pos);
    }
}

//---------------------------------------------------------------------------//

TEST_F(RectArrayProtoTest, errors)
{
    // Non-cuboid shape
    {
        PlacedShape::Params ps_params;
        ps_params.md     = ORANGE_MD_FROM_SOURCE("sphere");
        ps_params.shape  = std::make_shared<SphereShape>(3.0);
        auto sphere_unit = this->make_simple_unit(
            ORANGE_MD_FROM_SOURCE("sphere"), std::move(ps_params));

        RectArrayProto::Params params;
        params.units = ProtoVec3{{1, 1, 1}, sphere_unit};
        params.md    = ORANGE_MD_FROM_SOURCE("array of spheres???");
        EXPECT_THROW(RectArrayProto(std::move(params)), validation_error);
    }

    // Inconsistent widths
    {
        auto z3 = this->make_array_cell(ORANGE_MD_FROM_SOURCE("height=3"),
                                        {1, 2, 3});
        auto z1 = this->make_array_cell(ORANGE_MD_FROM_SOURCE("height=1"),
                                        {1, 2, 1});
        RectArrayProto::Params params;
        params.units = ProtoVec3{{1, 3, 1}, {z3, z1, z3}};
        params.md    = ORANGE_MD_FROM_SOURCE("bad widths");
        EXPECT_THROW(RectArrayProto(std::move(params)), validation_error);
    }

    // Missing important units
    {
        auto u
            = this->make_array_cell(ORANGE_MD_FROM_SOURCE("unit"), {1, 2, 3});
        RectArrayProto::Params params;
        params.units = ProtoVec3{{1, 1, 4}, {u, u, nullptr, u}};
        params.md    = ORANGE_MD_FROM_SOURCE("missing width");
        EXPECT_THROW(RectArrayProto(std::move(params)), validation_error);
    }

    // Rotated
    {
    }
}
