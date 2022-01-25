//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/tstUniverse.cc
 * \brief Tests for class Universe
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Universe.hh"

#include "celeritas_test.hh"
#include "base/FixedViewArray.hh"
#include "base/Future.hh"
#include "base/VectorFunctions.hh"
#include "orange/construct/ShapeContainer.hh"
#include "orange/construct/SphereShape.hh"
#include "orange/construct/UnitBuilder.hh"
#include "../SimpleUnitTracker.hh"

using namespace celeritas;
using make_fixed_view;
using std::make_shared;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class UniverseTest : public Test
{
  protected:
    //@{
    //! Type aliases
    using UPConstTracker  = std::unique_ptr<const celeritas::Tracker>;
    using SPConstUniverse = std::shared_ptr<const celeritas::Universe>;
    using State_t         = TrackingState;
    //@}

  protected:
    UPConstTracker build_spheres_tracker(real_type      radius,
                                         Real3          center,
                                         ObjectMetadata md) const
    {
        ShapeContainer shapes;

        // Build two spheres with scaled radius at the given position
        auto outer = shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("outer"), Transform{center}, radius);
        auto inner = shapes.emplace<SphereShape>(
            ORANGE_MD_FROM_SOURCE("inner"), Transform{center}, radius * 0.5);

        // Two interior regions
        UnitBuilder build;
        build.exterior(
            {{neg, outer}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("exterior"));
        build.region({{neg, outer}, {pos, inner}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("outer"));
        build.region(
            {{neg, inner}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("inner"));
        EXPECT_TRUE(build.is_simple());

        auto built = build(std::move(md));

        return make_unique<SimpleUnitTracker>(std::move(built.surfaces),
                                              std::move(built.regions));
    }
};

//---------------------------------------------------------------------------//
// SINGLE LEVEL
//---------------------------------------------------------------------------//
class SingleLevel : public UniverseTest
{
  protected:
    void SetUp()
    {
        Universe::Params params;
        params.id      = UniverseId{0};
        params.tracker = this->build_spheres_tracker(
            1.0, {0, 0, 0}, ORANGE_MD_FROM_SOURCE("single"));
        uptr = make_shared<Universe>(std::move(params));
    }

    SPConstUniverse uptr;
};

TEST_F(SingleLevel, accessors)
{
    const Universe& u = *uptr;
    EXPECT_EQ(UniverseId(0), u.id());
    EXPECT_EQ(3, u.num_volumes());
    EXPECT_EQ(2, u.num_surfaces());
    EXPECT_EQ(nullptr, u.daughter(VolumeId{0}));
}

TEST_F(SingleLevel, tracking)
{
    Real3                     pos{0, 0, 0};
    Real3                     dir{1, 0, 0};
    celeritas::SenseContainer temp_senses;
    celeritas::VecNextFace     temp_face_dist;

    LocalState state;
    state.pos            = make_fixed_view(pos);
    state.dir            = make_fixed_view(dir);
    state.cell           = {};
    state.surface        = {};
    state.sense          = {};
    state.temp_senses    = &temp_senses;
    state.temp_face_dist = &temp_face_dist;

    const Universe& u    = *uptr;
    auto            init = u.initialize(state);
    EXPECT_EQ(VolumeId{2}, init.cell);
    EXPECT_EQ(SurfaceId{}, init.surface);
    state.cell = init.cell;

    auto next = u.intersect(state);
    EXPECT_EQ(SurfaceId{1}, next.surface);
    EXPECT_SOFT_EQ(0.5, next.distance);
    state.surface = next.surface;
    state.sense   = next.sense;
    axpy(next.distance, state.dir, make_fixed_view(pos));

    Real3 norm = u.normal(state);
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), norm);
}
