//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstHexArrayProto.cc
 * \brief Tests for class HexArrayProto
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../HexArrayProto.hh"

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "base/Range.hh"
#include "base/StringFunctions.hh"
#include "orange/query/HexArrayMetadata.hh"
#include "orange/track/Definitions.hh"
#include "orange/track/HexArrayTracker.hh"
#include "../CuboidShape.hh"
#include "../PrismShape.hh"
#include "../PlacedShape.hh"
#include "../SphereShape.hh"
#include "../UnitProto.hh"
#include "ProtoTest.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class HexArrayProtoTest : public celeritas::ProtoTest
{
  protected:
    //// TYPE ALIASES ////
    using SPConstArray = std::shared_ptr<const ArrayProto>;
    using ProtoVec3    = HexArrayProto::ProtoVec3;
    using DimVector    = HexArrayProto::DimVector;
    using Layout       = HexArrayProto::Layout;
    using Orientation  = HexArrayProto::Orientation;
    using VecProto     = std::vector<SPConstProto>;

  public:
    // Make a simple unit with an untransformed hex prism
    SPConstProto make_array_cell(ObjectMetadata unit_md,
                                 real_type      apothem,
                                 Orientation    orientation,
                                 real_type      z_lower,
                                 real_type      z_upper) const
    {
        real_type rotation
            = (orientation == Orientation::pointy_top ? 0.5 : 0.0);
        PlacedShape::Params ps_params;
        ps_params.transform = {};
        ps_params.md        = ORANGE_MD_FROM_SOURCE("prism");
        ps_params.shape     = std::make_shared<PrismShape>(
            6, apothem, rotation, z_lower, z_upper);
        return this->make_simple_unit(std::move(unit_md), ps_params);
    }

    // Make cells with different universe names, all z from zero to 1
    VecProto
    make_several_cells(int count, real_type apothem, Orientation orientation)
    {
        CELER_EXPECT(count > 0);
        VecProto daughters(count);
        for (auto i : range(count))
        {
            daughters[i] = this->make_array_cell(
                ORANGE_MD_FROM_SOURCE("u" + std::to_string(i)),
                apothem,
                orientation,
                0.0,
                1.0);
        }
        return daughters;
    }
};

struct UnpackedBuildResult
{
    std::shared_ptr<const HexArrayTracker>  tracker;
    std::shared_ptr<const HexArrayMetadata> md;

    std::vector<std::string> fill;        //!< Fill names
    std::vector<real_type>   translation; //!< Flattened cell translations
    std::vector<real_type>   bbox; //!< Flattened bounding box {lower, upper}

    void print_expected() const;
};

void UnpackedBuildResult::print_expected() const
{
    cout
        << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
        << "static const char* const expected_fill[] = "
        << to_string(this->fill) << ";\n"
        << "EXPECT_VEC_EQ(expected_fill, built.fill);\n"
        << "static const real_type expected_translation[] = "
        << to_string(this->translation) << ";\n"
        << "EXPECT_VEC_SOFT_EQ(expected_translation, built.translation);\n"
        << "static const real_type expected_bbox[] = " << to_string(this->bbox)
        << ";\n"
        << "EXPECT_VEC_SOFT_EQ(expected_bbox, built.bbox);\n"
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
            = std::dynamic_pointer_cast<const HexArrayTracker>(tracker);
        EXPECT_TRUE(result.tracker);
    }
    {
        // Save metadata
        result.md = std::dynamic_pointer_cast<const HexArrayMetadata>(build.md);
        EXPECT_TRUE(result.md);
    }
    {
        // Fills
        CELER_ASSERT(build.md);
        EXPECT_EQ(0, build.matids.size());
        for (const auto& cell_daughter : build.daughters)
        {
            unsigned int cell = cell_daughter.first.unchecked_get();
            CELER_ASSERT(cell_daughter.second.proto);
            result.fill.push_back(
                std::to_string(cell) + ":"
                + cell_daughter.second.proto->metadata().name());

            const Transform& transform = cell_daughter.second.transform;
            EXPECT_FALSE(transform.has_rotation());
            const Real3& translation = transform.translation();
            result.translation.insert(result.translation.end(),
                                      translation.begin(),
                                      translation.end());
        }
    }
    {
        // Bbox
        const auto& bbox = build.md->bbox();
        result.bbox.insert(
            result.bbox.end(), bbox.lower().begin(), bbox.lower().end());
        result.bbox.insert(
            result.bbox.end(), bbox.upper().begin(), bbox.upper().end());
    }

    return result;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(HexArrayProtoTest, single_flat)
{
    // Hex-prism specification
    real_type apothem = 2.0;
    real_type z_lower = 0.0;
    real_type z_upper = 2.0;

    Layout       layout      = Layout::rhomboidal;
    Orientation  orientation = Orientation::flat_top;
    SPConstArray arr;

    {
        auto daughter = this->make_array_cell(ORANGE_MD_FROM_SOURCE("unit"),
                                              apothem,
                                              Orientation::flat_top,
                                              z_lower,
                                              z_upper);
        HexArrayProto::Params params;
        params.units  = ProtoVec3{{1, 1, 1}, ProtoVec3::Storage_t{daughter}};
        params.md     = ORANGE_MD_FROM_SOURCE("parent array");
        params.layout = layout;
        params.orientation = orientation;

        arr = std::make_shared<HexArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    // Check proto metadata
    EXPECT_EQ("parent array", arr->metadata().name());

    // Check placement
    {
        // Transform is from array cell center, not daughter
        auto transform = arr->calc_placement({0, 0, 0}, Real3{1.73, 0, 0});
        EXPECT_FALSE(transform.has_rotation());
        EXPECT_VEC_SOFT_EQ(Real3(1.73, 0, 0), transform.translation());
    }

    // Build
    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr->build(args));

    static const char* const expected_fill[] = {"0:unit"};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    static const real_type expected_translation[] = {0, 0, 0};
    EXPECT_VEC_EQ(expected_translation, built.translation);
    static const real_type expected_bbox[]
        = {-2.309401076759, -2, 0, 2.3094010767585, 2, 2};
    EXPECT_VEC_SOFT_EQ(expected_bbox, built.bbox);

    // Check built metadata
    ASSERT_SP_TRUE(built.md);
    EXPECT_EQ("parent array", built.md->metadata().name());
}

//---------------------------------------------------------------------------//
// This runs a load of combinations primarily to make sure no DBC checks fail
TEST_F(HexArrayProtoTest, combinations)
{
    std::vector<real_type> lower_trans, upper_trans;
    for (auto layout : {Layout::rhomboidal, Layout::rectangular})
    {
        for (auto orient : {Orientation::pointy_top, Orientation::flat_top})
        {
            // Make a cell with an inner diameter and height of one
            real_type apothem = 0.5;
            auto      cell    = this->make_array_cell(
                ORANGE_MD_FROM_SOURCE("unit"), apothem, orient, 0.0, 1.0);
            for (auto nu : range(1, 6))
            {
                for (auto nv : range(1, 6))
                {
                    HexArrayProto::Params params;
                    params.units
                        = ProtoVec3{DimVector(nu, nv, 2),
                                    ProtoVec3::Storage_t(nu * nv * 2, cell)};
                    params.md          = ORANGE_MD_FROM_SOURCE("parent array");
                    params.layout      = layout;
                    params.orientation = orient;

                    // Create array
                    HexArrayProto    arr(std::move(params));
                    Proto::BuildArgs args;
                    args.implicit_boundary = true; // array is always truncated

                    // Transform lower left cell center to origin
                    Real3 t = arr.calc_placement(DimVector(0, 0, 0),
                                                 Real3(0, 0, 0))
                                  .translation();
                    if (layout == Layout::rhomboidal)
                    {
                        EXPECT_VEC_SOFT_EQ(Real3(0, 0, 0), t);
                    }

                    // Transform upper-right cell toward the upper right
                    t = arr.calc_placement(DimVector(nu - 1, nv - 1, 1),
                                           Real3(nu, nv, 1))
                            .translation();

                    // Build and check dimensions: number of built daughters
                    // should be equal to the number of user-input protos
                    auto built = unpack_proto_build(arr.build(args));
                    EXPECT_EQ(built.fill.size(), nu * nv * 2);
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
TEST_F(HexArrayProtoTest, z_translations)
{
    // Hex-prism specification
    real_type              apothem = 2.0;
    std::vector<real_type> z_edges = {0.0, 2.0, 5.0, 10.0};

    Layout       layout = Layout::rhomboidal;
    Orientation  orient = Orientation::flat_top;
    SPConstArray arr;

    {
        ProtoVec3::Storage_t daughters;
        for (auto i : range(3))
        {
            daughters.push_back(this->make_array_cell(
                ORANGE_MD_FROM_SOURCE("u" + std::to_string(i)),
                apothem,
                orient,
                z_edges[i],
                z_edges[i + 1]));
        }

        HexArrayProto::Params params;
        params.units  = ProtoVec3{{1, 1, 3}, daughters};
        params.md     = ORANGE_MD_FROM_SOURCE("parent array");
        params.layout = layout;

        arr = std::make_shared<HexArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    // Check placement
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
    static const char* const expected_fill[]      = {"0:u0", "1:u1", "2:u2"};
    static const real_type expected_translation[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_SOFT_EQ(expected_translation, built.translation);

    // Check built metadata
    ASSERT_SP_TRUE(built.md);
    EXPECT_EQ("parent array", built.md->metadata().name());
}

//---------------------------------------------------------------------------//
//     o   -
//
//   -   o
//
// -   o
TEST_F(HexArrayProtoTest, odd_pointy)
{
    HexArrayProto::Params params;
    params.units = ProtoVec3{
        {1, 3, 1}, this->make_several_cells(3, 1.0, Orientation::pointy_top)};
    params.md     = ORANGE_MD_FROM_SOURCE("parent array");
    params.layout = Layout::rectangular;

    HexArrayProto arr(std::move(params));
    EXPECT_EQ(DimVector(1, 0, 0), arr.user_to_tracker({0, 0, 0}));
    EXPECT_EQ(DimVector(1, 1, 0), arr.user_to_tracker({0, 1, 0}));
    EXPECT_EQ(DimVector(0, 2, 0), arr.user_to_tracker({0, 2, 0}));

    // Build
    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr.build(args));

    static const char* const expected_fill[] = {"3:u0", "4:u1", "2:u2"};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    static const real_type expected_translation[]
        = {2, 0, 0, 3, 1.732050807569, 0, 2, 3.464101615138, 0};
    EXPECT_VEC_SOFT_EQ(expected_translation, built.translation);
    static const real_type expected_bbox[]
        = {1, -1.154700538379, 0, 4, 4.618802153517, 1};
    EXPECT_VEC_SOFT_EQ(expected_bbox, built.bbox);
}

//---------------------------------------------------------------------------//
//       -
//    o
// o     o
//    -
// -
TEST_F(HexArrayProtoTest, odd_flat)
{
    HexArrayProto::Params params;
    params.units = ProtoVec3{
        {3, 1, 1}, this->make_several_cells(3, 1.0, Orientation::flat_top)};
    params.md     = ORANGE_MD_FROM_SOURCE("array");
    params.layout = Layout::rectangular;

    HexArrayProto arr(std::move(params));
    EXPECT_EQ(DimVector(0, 1, 0), arr.user_to_tracker({0, 0, 0}));
    EXPECT_EQ(DimVector(1, 1, 0), arr.user_to_tracker({1, 0, 0}));
    EXPECT_EQ(DimVector(2, 0, 0), arr.user_to_tracker({2, 0, 0}));

    // Build
    Proto::BuildArgs args;
    args.implicit_boundary = true; // array is always truncated
    auto built             = unpack_proto_build(arr.build(args));

    static const char* const expected_fill[] = {"1:u0", "3:u1", "4:u2"};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    static const real_type expected_translation[]
        = {0, 2, 0, 1.732050807569, 3, 0, 3.464101615138, 2, 0};
    EXPECT_VEC_SOFT_EQ(expected_translation, built.translation);
    static const real_type expected_bbox[]
        = {-1.154700538379, 1, 0, 4.618802153517, 4, 1};
    EXPECT_VEC_SOFT_EQ(expected_bbox, built.bbox);
}

//---------------------------------------------------------------------------//
TEST_F(HexArrayProtoTest, daughter_layout)
{
    SPConstArray arr;
    {
        HexArrayProto::Params params;
        params.units = ProtoVec3{
            {4, 2, 1},
            this->make_several_cells(4 * 2, 2.0, Orientation::flat_top)};
        params.md     = ORANGE_MD_FROM_SOURCE("parent array");
        params.layout = Layout::rhomboidal;
        arr           = std::make_shared<HexArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    // Check placement
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
    static const char* const expected_fill[]
        = {"0:u0", "1:u1", "2:u2", "3:u3", "4:u4", "5:u5", "6:u6", "7:u7"};
    // clang-format off
    static const real_type expected_translation[] = {
        0.0,            0.0,  0.0,
        0.0,            4.0,  0.0,
        3.464101615138, 2.0,  0.0,
        3.464101615138, 6.0,  0.0,
        6.928203230276, 4.0,  0.0,
        6.928203230276, 8.0,  0.0,
        10.39230484541, 6.0,  0.0,
        10.39230484541, 10.0, 0.0};
    // clang-format on
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_SOFT_EQ(expected_translation, built.translation);
    static const real_type expected_bbox[]
        = {-2.309401076759, -2., 0, 12.7017059221718, 12, 1};
    EXPECT_VEC_SOFT_EQ(expected_bbox, built.bbox);

    // Check built metadata
    ASSERT_SP_TRUE(built.md);
    EXPECT_EQ("parent array", built.md->metadata().name());
}

//---------------------------------------------------------------------------//

TEST_F(HexArrayProtoTest, errors)
{
    // Hex-prism specification
    real_type      apothem = 2.0;
    real_type      z_lower = 0.0;
    real_type      z_upper = 2.0;
    constexpr auto flat    = Orientation::flat_top;
    constexpr auto pointy  = Orientation::pointy_top;

    // Empty array
    {
        HexArrayProto::Params params;
        params.md     = ORANGE_MD_FROM_SOURCE("sphere_array");
        params.layout = Layout::rhomboidal;

        EXPECT_THROW(HexArrayProto(std::move(params)), assertion);
    }

    // Non-hex shape
    {
        PlacedShape::Params ps_params;
        ps_params.md     = ORANGE_MD_FROM_SOURCE("sphere");
        ps_params.shape  = std::make_shared<SphereShape>(3.0);
        auto sphere_unit = this->make_simple_unit(
            ORANGE_MD_FROM_SOURCE("sphere"), std::move(ps_params));

        HexArrayProto::Params params;
        params.units  = ProtoVec3{{1, 1, 1}, sphere_unit};
        params.md     = ORANGE_MD_FROM_SOURCE("sphere_array");
        params.layout = Layout::rhomboidal;

        EXPECT_THROW(HexArrayProto(std::move(params)), validation_error);
    }

    // Prism shape other than hexagon
    {
        PlacedShape::Params ps_params;
        ps_params.md = ORANGE_MD_FROM_SOURCE("prism");
        ps_params.shape
            = std::make_shared<PrismShape>(7, apothem, 0.0, z_lower, z_upper);
        auto prism_unit = this->make_simple_unit(
            ORANGE_MD_FROM_SOURCE("prism"), std::move(ps_params));

        HexArrayProto::Params params;
        params.units  = ProtoVec3{{1, 1, 1}, prism_unit};
        params.md     = ORANGE_MD_FROM_SOURCE("non_hex_prism_array");
        params.layout = Layout::rhomboidal;

        EXPECT_THROW(HexArrayProto(std::move(params)), validation_error);
    }

    // Prism shape (hexprism with rotation other than rotation = 0.0 and 0.5)
    {
        real_type rotation = 0.7;

        PlacedShape::Params ps_params;
        ps_params.md    = ORANGE_MD_FROM_SOURCE("prism");
        ps_params.shape = std::make_shared<PrismShape>(
            6, apothem, rotation, z_lower, z_upper);
        auto prism_unit = this->make_simple_unit(
            ORANGE_MD_FROM_SOURCE("prism"), std::move(ps_params));

        HexArrayProto::Params params;
        params.units = ProtoVec3{{1, 1, 1}, prism_unit};
        params.md = ORANGE_MD_FROM_SOURCE("hex_prism_with_invalid_rotation");
        params.layout = Layout::rhomboidal;

        EXPECT_THROW(HexArrayProto(std::move(params)), validation_error);
    }

    // Try to place flat_top hexprism to POINTY_TOP array
    {
        // Construct a unit with FLAT_TOP hexprism
        auto ft = this->make_array_cell(
            ORANGE_MD_FROM_SOURCE("flat_top"), 1.0, flat, z_lower, z_upper);

        HexArrayProto::Params params;
        params.units = ProtoVec3{{1, 1, 1}, ft};
        params.md
            = ORANGE_MD_FROM_SOURCE("place_pointy_top_into_flat_top_array");
        params.layout      = Layout::rhomboidal;
        params.orientation = pointy;

        EXPECT_THROW(HexArrayProto(std::move(params)), validation_error);
    }

    // Try to place both FLAT_TOP and POINTY_TOP hexprisms into an array whose
    // orientation has not been specified (In SCALE implementation, orientation
    // of array is always defined. However, omnibus side has this flexibility)
    {
        // Construct a unit with POINTY_TOP hexprism
        auto pt = this->make_array_cell(
            ORANGE_MD_FROM_SOURCE("pointy_top"), 1.0, pointy, z_lower, z_upper);

        // Construct a unit with FLAT_TOP hexprism
        auto ft = this->make_array_cell(
            ORANGE_MD_FROM_SOURCE("flat_top"), 1.0, flat, z_lower, z_upper);

        HexArrayProto::Params params;
        params.units  = ProtoVec3{{1, 2, 1}, {pt, ft}};
        params.md     = ORANGE_MD_FROM_SOURCE("bad orientation");
        params.layout = Layout::rhomboidal;

        EXPECT_THROW(HexArrayProto(std::move(params)), validation_error);
    }

    // Inconsistent apothems
    {
        auto a3 = this->make_array_cell(
            ORANGE_MD_FROM_SOURCE("apothem=3"), 3.0, flat, z_lower, z_upper);

        auto a1 = this->make_array_cell(
            ORANGE_MD_FROM_SOURCE("apothem=1"), 1.0, flat, z_lower, z_upper);

        HexArrayProto::Params params;
        params.units  = ProtoVec3{{1, 2, 1}, {a1, a3}};
        params.md     = ORANGE_MD_FROM_SOURCE("bad apothems");
        params.layout = Layout::rhomboidal;

        EXPECT_THROW(HexArrayProto(std::move(params)), validation_error);
    }

    // Inconsistent z-thickness
    {
        auto a3 = this->make_array_cell(
            ORANGE_MD_FROM_SOURCE("height=2"), 1.0, flat, z_lower, z_upper);

        auto a1 = this->make_array_cell(
            ORANGE_MD_FROM_SOURCE("height=1"), 1.0, flat, 1.0, z_upper);

        HexArrayProto::Params params;
        params.units  = ProtoVec3{{1, 2, 1}, {a1, a3}};
        params.md     = ORANGE_MD_FROM_SOURCE("bad_heights");
        params.layout = Layout::rhomboidal;

        EXPECT_THROW(HexArrayProto(std::move(params)), validation_error);
    }

    // Missing important units
    {
        auto u = this->make_array_cell(
            ORANGE_MD_FROM_SOURCE("unit"), 2, flat, z_lower, z_upper);
        HexArrayProto::Params params;
        params.units  = ProtoVec3{{1, 1, 4}, {u, u, nullptr, u}};
        params.md     = ORANGE_MD_FROM_SOURCE("missing_unit");
        params.layout = Layout::rhomboidal;

        EXPECT_THROW(HexArrayProto(std::move(params)), validation_error);
    }
}
