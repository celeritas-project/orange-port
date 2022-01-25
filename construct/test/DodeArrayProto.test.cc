//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstDodeArrayProto.cc
 * \brief Tests for class DodeArrayProto
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../DodeArrayProto.hh"

#include <memory>

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "base/StringFunctions.hh"
#include "orange/query/DodeArrayMetadata.hh"
#include "orange/track/DodeArrayTracker.hh"
#include "../CuboidShape.hh"
#include "../RhombicDodecahedronShape.hh"
#include "../PlacedShape.hh"
#include "../SphereShape.hh"
#include "../UnitProto.hh"
#include "ProtoTest.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class DodeArrayProtoTest : public celeritas::ProtoTest
{
  protected:
    //// TYPE ALIASES ////
    using SPConstProto = std::shared_ptr<const Proto>;
    using SPConstArray = std::shared_ptr<const ArrayProto>;
    using ProtoVec3    = DodeArrayProto::ProtoVec3;

  public:
    real_type                  apothem_ = 1.0;
    constexpr static real_type sqrt2    = constants::sqrt_two;

    // Make a simple unit with an untransformed dodeca
    SPConstProto
    make_array_cell(ObjectMetadata unit_md, real_type apothem) const
    {
        PlacedShape::Params ps_params;
        ps_params.transform = {};
        ps_params.md        = ORANGE_MD_FROM_SOURCE("dode");
        ps_params.shape = std::make_shared<RhombicDodecahedronShape>(apothem);
        return this->make_simple_unit(std::move(unit_md), ps_params);
    }
};

struct UnpackedBuildResult
{
    std::shared_ptr<const DodeArrayTracker>  tracker;
    std::shared_ptr<const DodeArrayMetadata> md;

    std::vector<std::string> fill;        //!< Fill names
    std::vector<real_type>   translation; //!< Flattened cell translations

    void print_expected() const;
};

void UnpackedBuildResult::print_expected() const
{
    cout << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
         << "static const char* const expected_fill[] = "
         << to_string(this->fill) << ";\n"
         << "static const real_type expected_translation[] = "
         << to_string(this->translation) << ";\n"
         << "EXPECT_VEC_EQ(expected_fill, built.fill);\n"
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
            = std::dynamic_pointer_cast<const DodeArrayTracker>(tracker);
        EXPECT_TRUE(result.tracker);
    }
    {
        // Save metadata
        result.md
            = std::dynamic_pointer_cast<const DodeArrayMetadata>(build.md);
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

    return result;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(DodeArrayProtoTest, single_unit)
{
    this->apothem_ = 1.0;
    SPConstArray arr;
    {
        auto daughter = this->make_array_cell(ORANGE_MD_FROM_SOURCE("unit"),
                                              this->apothem_);
        DodeArrayProto::Params params;
        params.units = ProtoVec3{{1, 1, 1}, ProtoVec3::Storage_t{daughter}};
        params.md    = ORANGE_MD_FROM_SOURCE("parent array");

        arr = std::make_shared<DodeArrayProto>(std::move(params));
    }
    CELER_ASSERT(arr);

    // Check proto metadata
    EXPECT_EQ("parent array", arr->metadata().name());

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
    static const real_type   expected_translation[] = {0, 0, 0};
    EXPECT_VEC_EQ(expected_fill, built.fill);
    EXPECT_VEC_EQ(expected_translation, built.translation);

    // Check built metadata
    ASSERT_SP_TRUE(built.md);
    EXPECT_EQ("parent array", built.md->metadata().name());
}

//---------------------------------------------------------------------------//

TEST_F(DodeArrayProtoTest, errors)
{
    // Empty array
    {
        DodeArrayProto::Params params;
        params.md = ORANGE_MD_FROM_SOURCE("sphere_array");
        EXPECT_THROW(DodeArrayProto(std::move(params)), assertion);
    }

    // Non-dodeca shape
    {
        PlacedShape::Params ps_params;
        ps_params.md     = ORANGE_MD_FROM_SOURCE("sphere");
        ps_params.shape  = std::make_shared<SphereShape>(3.0);
        auto sphere_unit = this->make_simple_unit(
            ORANGE_MD_FROM_SOURCE("sphere"), std::move(ps_params));

        DodeArrayProto::Params params;
        params.units = ProtoVec3{{1, 1, 1}, sphere_unit};
        params.md    = ORANGE_MD_FROM_SOURCE("sphere_array");
        EXPECT_THROW(DodeArrayProto(std::move(params)), validation_error);
    }

    // Inconsistent apothems
    {
        auto a3 = this->make_array_cell(ORANGE_MD_FROM_SOURCE("a3"), 3);
        auto a1 = this->make_array_cell(ORANGE_MD_FROM_SOURCE("a1"), 1);
        DodeArrayProto::Params params;
        params.units = ProtoVec3{{1, 2, 1}, {a1, a3}};
        params.md    = ORANGE_MD_FROM_SOURCE("bad widths");
        EXPECT_THROW(DodeArrayProto(std::move(params)), validation_error);
    }

    // Missing important units
    {
        auto u = this->make_array_cell(ORANGE_MD_FROM_SOURCE("unit"), 2);
        DodeArrayProto::Params params;
        params.units = ProtoVec3{{1, 1, 4}, {u, u, nullptr, u}};
        params.md    = ORANGE_MD_FROM_SOURCE("missing width");
        EXPECT_THROW(DodeArrayProto(std::move(params)), assertion);
    }
}
