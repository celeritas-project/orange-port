//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/tstTrackingGeometry.cc
 * \brief Tests for class TrackingGeometry
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../TrackingGeometry.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "base/FixedViewArray.hh"
#include "base/Future.hh"
#include "base/VectorFunctions.hh"
#include "orange/TransformUtils.hh"
#include "orange/Fuzziness.hh"
#include "orange/construct/CylinderShape.hh"
#include "orange/construct/IntersectionShape.hh"
#include "orange/construct/PlaneShape.hh"
#include "orange/construct/ShapeContainer.hh"
#include "orange/construct/SphereShape.hh"
#include "orange/construct/UnitBuilder.hh"
#include "orange/query/UnitMetadata.hh"
#include "../SimpleUnitTracker.hh"
#include "../StateView.hh"
#include "../Universe.hh"

using namespace celeritas;

//---------------------------------------------------------------------------//
// Test fixture base class
class TrackingGeometryTest : public ::Test
{
  protected:
    //@{
    //! Type aliases
    using SPConstUniverse    = std::shared_ptr<const Universe>;
    using UPConstTracker     = std::unique_ptr<const celeritas::Tracker>;
    using SPConstUnit_MD     = std::shared_ptr<const UnitMetadata>;
    using Daughter           = celeritas::Universe::Daughter;
    using MapSurfaceBoundary = TrackingGeometry::MapSurfaceBoundary;
    //@}

    struct TrackResult
    {
        std::vector<std::string> cells;
        std::vector<std::string> surfaces;
        std::vector<char>        senses;
        std::vector<real_type>   distances;
        std::vector<real_type>   normals;

        void print_expected() const;
    };

    static constexpr real_type inf = HUGE_VAL;

  protected:
    //// CONSTRUCTION ////

    void SetUp() override
    {
        celeritas::fuzziness() = Fuzziness{};

        this->build_universes();

        CELER_EXPECT(!universes_.empty());
        CELER_EXPECT(std::all_of(
            universes_.begin(), universes_.end(), [](const SPConstUniverse& u) {
                return static_cast<bool>(u);
            }));

        TrackingGeometry::Params params;
        params.universes.assign(universes_.begin(), universes_.end());
        params.boundaries = this->build_boundaries();

        geo_ = std::make_shared<TrackingGeometry>(std::move(params));
        CELER_ENSURE(geo_);
    }

    // Construct the geometry (daughter must implement)
    virtual void build_universes() = 0;

    // Create boundaries
    virtual MapSurfaceBoundary build_boundaries() const { return {}; }

    // Build two concentric spheres, with an optional daughter for the inner
    // sphere (with half the radius of the parent).
    SPConstUniverse build_spheres_universe(real_type      radius,
                                           Real3          center,
                                           Daughter       daughter,
                                           UniverseId     id,
                                           ObjectMetadata md);

    // Build a single sphere at the center.
    SPConstUniverse build_sphere_univ(real_type      radius,
                                      Daughter       daughter,
                                      UniverseId     id,
                                      ObjectMetadata md);

    void insert_universe(SPConstUniverse u, SPConstUnit_MD md);

    //// ACCESSORS ////

    // Get the built geometry
    const TrackingGeometry& geo() const { return *geo_; }

    // Print all metadata
    void describe(std::ostream& os) const;

    //// TESTING ////

    // Track until the geometry is exited
    TrackResult track(SpanConstReal3 pos, Real3 dir);

  private:
    //// DATA ////
    std::vector<SPConstUniverse>      universes_;
    std::vector<SPConstUnit_MD>       md_;
    std::shared_ptr<TrackingGeometry> geo_;
};

constexpr real_type TrackingGeometryTest::inf;

//---------------------------------------------------------------------------//
/*!
 * Build two concentric spheres.
 */
auto TrackingGeometryTest::build_spheres_universe(real_type      radius,
                                                  Real3          center,
                                                  Daughter       daughter,
                                                  UniverseId     id,
                                                  ObjectMetadata md)
    -> SPConstUniverse
{
    CELER_EXPECT(radius > 0);
    ShapeContainer shapes;

    // Build two spheres with scaled radius at the given position
    auto outer = shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("outer"), Transform{center}, radius);
    auto inner = shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("inner"), Transform{center}, radius * 0.5);

    auto exterior_zorder
        = (id == UniverseId{0} ? ZOrder::media : ZOrder::implicit_exterior);

    // Two interior regions
    UnitBuilder build;
    build.exterior(
        {{neg, outer}}, exterior_zorder, ORANGE_MD_FROM_SOURCE("exterior"));
    build.region({{neg, outer}, {pos, inner}},
                 ZOrder::media,
                 ORANGE_MD_FROM_SOURCE("outer"));
    build.region({{neg, inner}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("inner"));

    // Construct regions
    auto built = build(std::move(md));
    EXPECT_TRUE(built.md->is_simple());

    // Construct universe params
    Universe::Params params;
    params.id      = id;
    params.tracker = make_unique<SimpleUnitTracker>(std::move(built.surfaces),
                                                    std::move(built.regions));
    if (daughter.universe)
    {
        // Inner region is the daughter
        params.daughters.resize(params.tracker->num_volumes());
        params.daughters.back() = std::move(daughter);
    }

    auto u = std::make_shared<Universe>(std::move(params));
    this->insert_universe(u, std::move(built.md));
    return u;
}

//---------------------------------------------------------------------------//
/*!
 * Build a single sphere.
 */
auto TrackingGeometryTest::build_sphere_univ(real_type      radius,
                                             Daughter       daughter,
                                             UniverseId     id,
                                             ObjectMetadata md)
    -> SPConstUniverse
{
    auto sph = std::make_shared<PlacedShape>(
        PlacedShape::Params{std::make_shared<SphereShape>(radius),
                            {},
                            ORANGE_MD_FROM_SOURCE("sph")});

    auto exterior_zorder
        = (id == UniverseId{0} ? ZOrder::media : ZOrder::implicit_exterior);

    UnitBuilder build;
    build.exterior(
        {{neg, sph}}, exterior_zorder, ORANGE_MD_FROM_SOURCE("exterior"));
    build.region(
        {{neg, sph}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("interior"));

    // Construct regions
    auto built = build(std::move(md));
    EXPECT_TRUE(built.md->is_simple());

    // Construct universe params
    Universe::Params params;
    params.id      = id;
    params.tracker = make_unique<SimpleUnitTracker>(std::move(built.surfaces),
                                                    std::move(built.regions));
    if (daughter.universe)
    {
        // Inner region is the daughter
        params.daughters.resize(params.tracker->num_volumes());
        params.daughters.back() = std::move(daughter);
    }

    auto u = std::make_shared<Universe>(std::move(params));
    this->insert_universe(u, std::move(built.md));
    return u;
}

void TrackingGeometryTest::insert_universe(SPConstUniverse u, SPConstUnit_MD md)
{
    CELER_EXPECT(u);
    auto id = u->id();
    CELER_EXPECT(id);
    CELER_EXPECT(id.get() >= universes_.size() || !universes_[id.get()]);
    CELER_EXPECT(md);
    // Update capacity
    if (id.get() >= universes_.size())
    {
        universes_.resize(id.get() + 1);
        md_.resize(id.get() + 1);
    }

    md_[id.get()]        = md;
    universes_[id.get()] = u;

    CELER_ENSURE(universes_[id.get()]);
    CELER_ENSURE(md_[id.get()]);
}

//---------------------------------------------------------------------------//
/*!
 * Print metadata to screen.
 */
void TrackingGeometryTest::describe(std::ostream& os) const
{
    auto get_daughter_md = [this](UniverseId uid) -> const ObjectMetadata& {
        CELER_EXPECT(uid < md_.size());
        CELER_ENSURE(md_[uid.get()]);
        return md_[uid.get()]->metadata();
    };
    auto print_local_vol = [](std::ostream& o, VolumeId i) {
        o << "(local cell " << i.get() << ")";
    };

    for (auto i : range(universes_.size()))
    {
        os << "========================================\n"
           << "Universe " << i << ": " << md_[i]->metadata()
           << "\n"
              "========================================\n"
           << to_stream(
                  *md_[i], *universes_[i], get_daughter_md, print_local_vol);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Track until the geometry is exited.
 */
auto TrackingGeometryTest::track(SpanConstReal3 pos, Real3 dir) -> TrackResult
{
    CELER_EXPECT(geo_);

    TrackingState   state;
    const StateView state_view(state);
    TrackResult     result;

    bool failed = false;
    normalize_direction(&dir);
    try
    {
        geo_->initialize(state, pos, dir);
    }
    catch (const std::exception& e)
    {
        ADD_FAILURE() << "Error during initialize: " << e.what();
        result.cells.push_back("[error]");
        return result;
    }

    // If outside, initial ray trace to inside of geometry
    if (geo_->boundary_state(state) != geometria::inside)
    {
        result.cells.push_back("---");

        // Now transport to the inside if possible
        real_type distance = geo_->distance_to_boundary(state);
        result.distances.push_back(distance);

        if (!std::isinf(distance))
        {
            // Finite distance to interior boundary
            try
            {
                geo_->cross_surface(state);
            }
            catch (const std::exception& e)
            {
                ADD_FAILURE() << "Error during cross_surface: " << e.what();
                failed = true;
            }
        }
    }

    while (geo_->boundary_state(state) == geometria::inside)
    {
        CELER_ASSERT(state.level >= 0);
        const UnitMetadata& md = *md_[state_view.universe().get()];

        result.cells.push_back(md.metadata().name() + ":"
                               + md.id_to_label(state_view.volume()));

        // Save beginning-of-segment data
        Real3       normal{0, 0, 0};
        std::string id_to_label("---");
        char        sense = ' ';
        if (auto univ_surf = state_view.surface())
        {
            const UnitMetadata& md = *md_[univ_surf.univ_id.get()];
            id_to_label            = md.id_to_label(state.on_surface.surface);
            sense                  = to_char(state.on_surface.sense);
            normal                 = geo_->normal(state);
        }
        result.surfaces.push_back(id_to_label);
        result.senses.push_back(sense);
        result.normals.insert(
            result.normals.end(), normal.begin(), normal.end());

        // Calculate next distance
        try
        {
            geo_->intersect(state);
        }
        catch (const std::exception& e)
        {
            ADD_FAILURE() << "Error during intersect: " << e.what();
            result.cells.push_back("[error]");
            failed = true;
            break;
        }
        result.distances.push_back(state_view.next_distance());

        // Move across surface
        try
        {
            geo_->cross_surface(state);
        }
        catch (const std::exception& e)
        {
            ADD_FAILURE() << "Error during cross_surface: " << e.what();
            result.cells.push_back("[error]");
            failed = true;
            break;
        }
    }

    if (!failed)
    {
        // Save end-of-segment surface data
        Real3       normal{0, 0, 0};
        std::string id_to_label("---");
        char        sense = ' ';
        if (auto univ_surf = state_view.surface())
        {
            const UnitMetadata& md = *md_[univ_surf.univ_id.get()];
            id_to_label            = md.id_to_label(state.on_surface.surface);
            sense                  = to_char(state.on_surface.sense);
            normal                 = geo_->normal(state);
        }
        result.senses.push_back(sense);
        result.surfaces.push_back(id_to_label);
        result.normals.insert(
            result.normals.end(), normal.begin(), normal.end());
    }

    return result;
}

void TrackingGeometryTest::TrackResult::print_expected() const
{
    cout << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
         << "static const char* const expected_cells[] = "
         << to_string(this->cells) << ";\n"
         << "EXPECT_VEC_EQ(expected_cells, result.cells);\n"
         << "static const char* const expected_surfaces[] = "
         << to_string(this->surfaces) << ";\n"
         << "EXPECT_VEC_EQ(expected_surfaces, result.surfaces);\n"
         << "static const char expected_senses[] = " << to_string(this->senses)
         << ";\n"
         << "EXPECT_VEC_EQ(expected_senses, result.senses);\n"
         << "static const real_type expected_distances[] = "
         << to_string(this->distances) << ";\n"
         << "EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);\n"
         << "static const real_type expected_normals[] = "
         << to_string(this->normals) << ";\n"
         << "EXPECT_VEC_SOFT_EQ(expected_normals, result.normals);\n"
         << "/*** END CODE ***/\n";
}

//---------------------------------------------------------------------------//
// Single-level case
class SingleLevelTest : public TrackingGeometryTest
{
    void build_universes() final
    {
        this->build_spheres_universe(
            1.0, {0, 0, 0}, {}, UniverseId{0}, ORANGE_MD_FROM_SOURCE("world"));
    }
};

//---------------------------------------------------------------------------//

TEST_F(SingleLevelTest, accessors)
{
    // this->describe(cout);
    EXPECT_EQ(1, geo().num_universes());
    EXPECT_EQ(UniverseId{0}, geo().get(UniverseId{0}).id());
}

//---------------------------------------------------------------------------//

TEST_F(SingleLevelTest, manual_track)
{
    TrackingState           state;
    StateView               state_view(state);
    const TrackingGeometry& g = this->geo();

    // Initialize
    g.initialize(state, {0.1, 0, 0}, {1, 0, 0});
    ASSERT_EQ(0, state.level);

    const auto& level_state = state.level_state[0];
    EXPECT_VEC_SOFT_EQ(Real3(0.1, 0, 0), level_state.pos);
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), level_state.dir);
    EXPECT_EQ(VolumeId{2}, level_state.cell);
    EXPECT_EQ(UniverseId{0}, level_state.univ_id);
    EXPECT_FALSE(state.on_surface);
    EXPECT_EQ(-1, state.next_level);
    EXPECT_SOFT_EQ(0.0, state.movement);

    // Test non-changing accesssors
    EXPECT_EQ(SurfaceId{}, state_view.surface().surface);
    EXPECT_VEC_SOFT_EQ(Real3(0.1, 0, 0), state_view.calc_position());
    EXPECT_EQ(geometria::inside, g.boundary_state(state));

    // Find next surface
    EXPECT_SOFT_EQ(0.4, g.distance_to_boundary(state));
    // Find next surface (cached)
    EXPECT_SOFT_EQ(0.4, g.distance_to_boundary(state));
    EXPECT_SOFT_EQ(0.4, state_view.next_distance());

    // Move inside
    state_view.move_internal(0.3);
    EXPECT_VEC_SOFT_EQ(Real3(0.4, 0, 0), state_view.calc_position());
    EXPECT_SOFT_EQ(0.1, g.distance_to_boundary(state));
    EXPECT_SOFT_EQ(0.1, state_view.next_distance());

    // Change direction
    g.set_direction(state, Real3{0, 1, 0});
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), level_state.dir);
    EXPECT_EQ(0.0, state.movement);
    EXPECT_VEC_SOFT_EQ(Real3(0.4, 0, 0), state_view.calc_position());
    // Find next surface
    g.intersect(state);
    EXPECT_SOFT_EQ(0.3, state_view.next_distance());

    // Next surface
    g.cross_surface(state);
    ASSERT_EQ(0, state.level);
    EXPECT_VEC_SOFT_EQ(Real3(0.4, 0.3, 0), level_state.pos);
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), level_state.dir);
    EXPECT_EQ(VolumeId{1}, level_state.cell);
    EXPECT_EQ(UniverseId{0}, level_state.univ_id);
    EXPECT_TRUE(state.on_surface);
    EXPECT_FALSE(state_view.has_next_distance());
    EXPECT_SOFT_EQ(0.0, state.movement);

    EXPECT_EQ(geometria::inside, g.boundary_state(state));
    EXPECT_EQ(SurfaceId{1}, state_view.surface().surface);

    // Move off surface
    EXPECT_SOFT_EQ(0.61651513899116805, g.distance_to_boundary(state));
    EXPECT_EQ(0, state.next_level);
    state_view.move_internal(0.1);
    EXPECT_EQ(SurfaceId{}, state_view.surface().surface);

    // Move outside
    g.cross_surface(state);
    EXPECT_EQ(geometria::outside, g.boundary_state(state));
}

//---------------------------------------------------------------------------//

TEST_F(SingleLevelTest, tracking)
{
    {
        // Outside heading in
        TrackResult result = this->track({-10, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"---", "world:outer", "world:inner", "world:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"outer.s", "inner.s", "inner.s", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {'-', '-', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {9, 0.5, 1, 0.5};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
        static const real_type expected_normals[]
            = {-1, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0};
        EXPECT_VEC_SOFT_EQ(expected_normals, result.normals);
    }

    {
        // Outside heading away
        TrackResult result = this->track({-10, 0, 0}, {-1, 0, 0});

        static const char* const expected_cells[] = {"---"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {"---"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' '};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Inside heading out
        TrackResult result = this->track({0, 0, 0}, {0, 1, 1});

        static const char* const expected_cells[]
            = {"world:inner", "world:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "inner.s", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {0.5, 0.5};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
        static const real_type expected_normals[] = {0,
                                                     0,
                                                     0,
                                                     0,
                                                     0.7071067811865,
                                                     0.7071067811865,
                                                     0,
                                                     0.7071067811865,
                                                     0.7071067811865};
        EXPECT_VEC_SOFT_EQ(expected_normals, result.normals);
    }

    {
        // Outside, not intersecting but inside bounding region
        TrackResult result = this->track({-9, 9, 0}, {1, 1, 0});

        static const char* const expected_cells[] = {"---"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const real_type expected_distances[] = {inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Outside, nowhere near
        TrackResult result = this->track({20, 20, 0}, {1, 0, 0});

        static const char* const expected_cells[] = {"---"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const real_type expected_distances[] = {inf};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

//---------------------------------------------------------------------------//
// Reflecting boundary
class ReflectingTest : public SingleLevelTest
{
    // Create boundaries
    MapSurfaceBoundary build_boundaries() const final
    {
        return {{SurfaceId{0}, geometria::REFLECT}};
    }
};

//---------------------------------------------------------------------------//

TEST_F(ReflectingTest, bounce_bounce_bounce)
{
    TrackingState           state;
    const TrackingGeometry& g = this->geo();

    g.initialize(state, {0.8, 0.2, 0}, {0, 1, 0});
    ASSERT_EQ(0, state.level);
    EXPECT_EQ(VolumeId{1}, state.level_state[0].cell);
    g.intersect(state);
    ASSERT_EQ(0, state.next_level);
    EXPECT_EQ(SurfaceId{0}, state.level_state[0].next.surface);
    g.cross_surface(state);
    EXPECT_EQ(VolumeId{0}, state.level_state[0].cell);
    EXPECT_EQ(SurfaceId{0}, state.on_surface.surface);
    EXPECT_EQ(geometria::REFLECT, g.boundary_state(state));
}

//---------------------------------------------------------------------------//
/*!
 * Multi-level coincident boundary test.
 *
 * Each level is rotated 90 degrees clockwise from the parent level, but all
 * levels are coincident.
 */

class MultiReflectingTest : public TrackingGeometryTest
{
    // Create boundaries
    MapSurfaceBoundary build_boundaries() const final
    {
        return {{SurfaceId{0}, geometria::REFLECT}};
    }

    void build_universes() final
    {
        Universe::Daughter daughter;

        daughter.universe = this->build_sphere_univ(
            1.0, daughter, UniverseId{2}, ORANGE_MD_FROM_SOURCE("princeps"));
        daughter.transform = geometria::rotation_matrix(Axis::z, -0.25);

        daughter.universe = this->build_sphere_univ(
            1.0, daughter, UniverseId{1}, ORANGE_MD_FROM_SOURCE("ochotona"));
        daughter.transform = geometria::rotation_matrix(Axis::z, -0.25);

        this->build_sphere_univ(
            1.0, daughter, UniverseId{0}, ORANGE_MD_FROM_SOURCE("lagomorpha"));
    }
};

//---------------------------------------------------------------------------//

TEST_F(MultiReflectingTest, bounce_bounce_bounce)
{
    TrackingState           state;
    const TrackingGeometry& g = this->geo();

    g.initialize(state, {0.8, 0.50, 0}, {0, 1, 0});
    ASSERT_EQ(2, state.level);
    EXPECT_EQ(VolumeId{1}, state.level_state[0].cell);
    EXPECT_EQ(VolumeId{1}, state.level_state[1].cell);
    EXPECT_EQ(VolumeId{1}, state.level_state[2].cell);
    EXPECT_VEC_SOFT_EQ(Real3(-0.8, -0.5, 0), state.level_state[2].pos);
    EXPECT_VEC_SOFT_EQ(Real3(0, -1, 0), state.level_state[2].dir);

    g.intersect(state);
    ASSERT_EQ(0, state.next_level);
    EXPECT_EQ(SurfaceId{0}, state.level_state[0].next.surface);
    EXPECT_SOFT_EQ(0.1, state.level_state[0].next.distance);

    g.cross_surface(state);
    ASSERT_EQ(0, state.level);
    EXPECT_EQ(VolumeId{0}, state.level_state[0].cell);
    EXPECT_EQ(SurfaceId{0}, state.on_surface.surface);
    EXPECT_EQ(geometria::REFLECT, g.boundary_state(state));
    EXPECT_VEC_SOFT_EQ(Real3(0.8, 0.6, 0), state.level_state[0].pos);

    g.reflect(state);
    EXPECT_EQ(2, state.level);
    EXPECT_EQ(VolumeId{1}, state.level_state[0].cell);
    EXPECT_EQ(VolumeId{1}, state.level_state[1].cell);
    EXPECT_EQ(VolumeId{1}, state.level_state[2].cell);
    EXPECT_VEC_SOFT_EQ(Real3(0.8, 0.6, 0), state.level_state[0].pos);
    EXPECT_VEC_SOFT_EQ(Real3(-0.6, 0.8, 0), state.level_state[1].pos);
    EXPECT_VEC_SOFT_EQ(Real3(-0.8, -0.6, 0), state.level_state[2].pos);
    EXPECT_VEC_SOFT_EQ(Real3(-.96, .28, 0), state.level_state[0].dir);
    EXPECT_VEC_SOFT_EQ(Real3(-.28, -.96, 0), state.level_state[1].dir);
    EXPECT_VEC_SOFT_EQ(Real3(.96, -.28, 0), state.level_state[2].dir);
}

//---------------------------------------------------------------------------//

TEST_F(MultiReflectingTest, internal_direction_change)
{
    TrackingState           state;
    StateView               state_view(state);
    const TrackingGeometry& g = this->geo();

    g.initialize(state, {0.8, 0.50, 0}, {0, 1, 0});
    g.intersect(state);
    ASSERT_EQ(2, state.level);
    EXPECT_SOFT_EQ(0.1, state.level_state[0].next.distance);

    state_view.move_internal(0.07);
    EXPECT_SOFT_EQ(0.03, g.distance_to_boundary(state));

    g.set_direction(state, {1, 0, 0});
    ASSERT_EQ(2, state.level);
    EXPECT_VEC_SOFT_EQ(Real3(0.8, 0.57, 0), state.level_state[0].pos);
    EXPECT_VEC_SOFT_EQ(Real3(-0.57, 0.8, 0), state.level_state[1].pos);
    EXPECT_VEC_SOFT_EQ(Real3(-0.8, -0.57, 0), state.level_state[2].pos);
    EXPECT_VEC_SOFT_EQ(Real3(1, 0, 0), state.level_state[0].dir);
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), state.level_state[1].dir);
    EXPECT_VEC_SOFT_EQ(Real3(-1, 0, 0), state.level_state[2].dir);

    EXPECT_SOFT_EQ(0.021644692065858728, g.distance_to_boundary(state));
}

//---------------------------------------------------------------------------//
// Three-level case
class ThreeLevelTest : public TrackingGeometryTest
{
    void build_universes() final
    {
        Universe::Daughter daughter;
        ObjectMetadata     md;

        // Inner
        md        = ORANGE_MD_FROM_SOURCE("babybear");
        auto baby = this->build_spheres_universe(
            0.5, {0.25, .5, 0}, daughter, UniverseId{2}, md);

        // Middle
        daughter.universe  = baby;
        daughter.transform = Transform{Real3{-.25, 0.5, 0}};
        md                 = ORANGE_MD_FROM_SOURCE("mamabear");
        auto mama          = this->build_spheres_universe(
            1.0, {0, 0.5, 0}, daughter, UniverseId{1}, md);

        // Outer
        daughter.universe  = mama;
        daughter.transform = Transform{Real3{0, -0.5, 0}};
        md                 = ORANGE_MD_FROM_SOURCE("papabear");
        this->build_spheres_universe(
            2.0, {0, 0.0, 0}, daughter, UniverseId{0}, md);
    }
};

//---------------------------------------------------------------------------//

TEST_F(ThreeLevelTest, internal_direction_change)
{
    TrackingState           state;
    StateView               state_view(state);
    const TrackingGeometry& g = this->geo();

    g.initialize(state, {0.1, 0.0, 0}, {1, 0, 0});
    ASSERT_EQ(2, state.level);
    EXPECT_VEC_SOFT_EQ(Real3(0.1, 0.0, 0), state.level_state[0].pos);
    EXPECT_VEC_SOFT_EQ(Real3(0.1, 0.5, 0), state.level_state[1].pos);
    EXPECT_VEC_SOFT_EQ(Real3(0.35, 0.0, 0), state.level_state[2].pos);

    g.intersect(state);
    EXPECT_EQ(1, state.next_level);
    EXPECT_SOFT_EQ(0.4, state.level_state[1].next.distance);

    state_view.move_internal(0.2);
    EXPECT_SOFT_EQ(0.2, g.distance_to_boundary(state));

    g.set_direction(state, {0, 1, 0});
    EXPECT_VEC_SOFT_EQ(Real3(0.3, 0.0, 0), state.level_state[0].pos);
    EXPECT_VEC_SOFT_EQ(Real3(0.3, 0.5, 0), state.level_state[1].pos);
    EXPECT_VEC_SOFT_EQ(Real3(0.55, 0.0, 0), state.level_state[2].pos);
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), state.level_state[0].dir);
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), state.level_state[1].dir);
    EXPECT_VEC_SOFT_EQ(Real3(0, 1, 0), state.level_state[2].dir);

    EXPECT_SOFT_EQ(0.4, g.distance_to_boundary(state));
}

TEST_F(ThreeLevelTest, tracking)
{
    {
        // Outside heading in
        TrackResult result = this->track({-3, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[] = {"---",
                                                     "papabear:outer",
                                                     "mamabear:outer",
                                                     "babybear:outer",
                                                     "mamabear:outer",
                                                     "papabear:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {
            "outer.s", "inner.s", "inner.s", "inner.s", "inner.s", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {'-', '-', '-', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1, 0.5, 1, 0.5, 1};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
        static const real_type expected_normals[]
            = {-1, 0, 0, -1, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0};
        EXPECT_VEC_SOFT_EQ(expected_normals, result.normals);
    }
    {
        // Inside heading out
        TrackResult              result = this->track({0, 0, 0}, {1, 0, 0});
        static const char* const expected_cells[]
            = {"babybear:outer", "mamabear:outer", "papabear:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "inner.s", "inner.s", "outer.s"};
        static const char expected_senses[] = {' ', '+', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const real_type expected_distances[] = {0.5, 0.5, 1};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
        static const real_type expected_normals[]
            = {0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0};
        EXPECT_VEC_SOFT_EQ(expected_normals, result.normals);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Nearly coincident boundary test.
 *
 * Each level is very slightly smaller than the parent level. This is a
 * surrogate for floating point errors that might crop up in array and hole
 * placement. Initialization should push the track to the interior surfaces
 */
class NearlyCoincidentTest : public TrackingGeometryTest
{
    // Create boundaries
    MapSurfaceBoundary build_boundaries() const final
    {
        return {{SurfaceId{0}, geometria::REFLECT}};
    }

    void build_universes() final
    {
        Universe::Daughter daughter;

        real_type eps = 1e-10;
        daughter.universe
            = this->build_sphere_univ(1.0 - 2 * eps,
                                      daughter,
                                      UniverseId{3},
                                      ORANGE_MD_FROM_SOURCE("glama"));

        daughter.universe = this->build_sphere_univ(
            1.0 - eps, daughter, UniverseId{2}, ORANGE_MD_FROM_SOURCE("lama"));

        daughter.universe = this->build_sphere_univ(
            1.0, daughter, UniverseId{1}, ORANGE_MD_FROM_SOURCE("camelidae"));

        this->build_spheres_universe(2.0,
                                     Real3{0, 0, 0},
                                     daughter,
                                     UniverseId{0},
                                     ORANGE_MD_FROM_SOURCE("mammalia"));
    }
};

//---------------------------------------------------------------------------//

TEST_F(NearlyCoincidentTest, tracking)
{
    {
        // Outside heading in
        TrackResult result = this->track({-10, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"---", "mammalia:outer", "glama:interior", "mammalia:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"outer.s", "inner.s", "inner.s", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {'-', '-', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {8, 1, 2, 1};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
        static const real_type expected_normals[]
            = {-1, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0};
        EXPECT_VEC_SOFT_EQ(expected_normals, result.normals);
    }
    {
        // Inside heading out
        TrackResult result = this->track({0, 0, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"glama:interior", "mammalia:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "inner.s", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 1};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
        static const real_type expected_normals[]
            = {0, 0, 0, -1, 0, 0, -1, 0, 0};
        EXPECT_VEC_SOFT_EQ(expected_normals, result.normals);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Skip-a-level coincident test.
 *
 * This is the case where the global unit has an array with shape A placed in a
 * hole with shape H; and inside array A is a single cell (for example) with an
 * object inside with shape H. The array correctly deletes its daughter cell
 * boundary, but can't delete the daughter shape H. So you end up starting on a
 * coincident surface.
 */
class GrandpaCoincidentTest : public TrackingGeometryTest
{
  protected:
    void build_universes() final;
};

void GrandpaCoincidentTest::build_universes()
{
    celeritas::fuzziness() = Fuzziness{1e-3};

    Universe::Daughter daughter;
    ShapeContainer     shapes;

    // Outer boundary
    auto outer = shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("outer"), Transform{}, 10.0);
    // Inner medium and hole shape
    auto inner = shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("inner"), Transform{}, 1.0);
    // Array shape (truncated by hole): cylinder with inf z extents
    auto array = shapes.emplace<CylinderShape>(
        ORANGE_MD_FROM_SOURCE("array"), Transform{}, 2.0, inf);

    // "Cell" unit: inner medium is same shape as placement hole in
    // grandparent.
    {
        UnitBuilder build;
        build.exterior({{inside, array}},
                       ZOrder::implicit_exterior,
                       ORANGE_MD_FROM_SOURCE("exterior"));
        build.region({{inside, array}, {outside, inner}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("midcell"));
        build.region(
            {{inside, inner}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("inner"));

        // Construct regions
        auto built = build(std::move(ORANGE_MD_FROM_SOURCE("baby")));
        EXPECT_TRUE(built.md->is_simple());

        // Construct universe params
        Universe::Params params;
        params.id      = UniverseId{2};
        params.tracker = make_unique<SimpleUnitTracker>(
            std::move(built.surfaces), std::move(built.regions));

        auto u            = std::make_shared<Universe>(std::move(params));
        daughter.universe = u;
        this->insert_universe(u, std::move(built.md));
    }

    // "Array" unit: in practice this would be a tracking array or
    // equivalent, but it *should* be representative for us to just make it a
    // larger shape than the parent hole. It's filled with a single cell that
    // has the same extents as it.
    {
        UnitBuilder build;
        build.exterior({{inside, array}},
                       ZOrder::implicit_exterior,
                       ORANGE_MD_FROM_SOURCE("exterior"));
        build.region({{inside, array}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("array cell 1"));

        // Construct regions
        auto built = build(std::move(ORANGE_MD_FROM_SOURCE("papa")));
        EXPECT_TRUE(built.md->is_simple());

        // Construct universe params
        Universe::Params params;
        params.id      = UniverseId{1};
        params.tracker = make_unique<SimpleUnitTracker>(
            std::move(built.surfaces), std::move(built.regions));
        // Inner region is the daughter
        params.daughters.resize(params.tracker->num_volumes());
        params.daughters.back() = std::move(daughter);

        auto u            = std::make_shared<Universe>(std::move(params));
        daughter.universe = u;
        this->insert_universe(u, std::move(built.md));
    }

    // Global unit
    {
        UnitBuilder build;
        build.exterior({{inside, outer}},
                       ZOrder::media,
                       ORANGE_MD_FROM_SOURCE("exterior"));
        build.region({{inside, outer}, {outside, inner}},
                     ZOrder::media,
                     ORANGE_MD_FROM_SOURCE("outer"));
        build.region(
            {{inside, inner}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("ahole"));

        // Construct regions
        auto built = build(std::move(ORANGE_MD_FROM_SOURCE("grandpa")));
        EXPECT_TRUE(built.md->is_simple());

        // Construct universe params
        Universe::Params params;
        params.id      = UniverseId{0};
        params.tracker = make_unique<SimpleUnitTracker>(
            std::move(built.surfaces), std::move(built.regions));
        // Inner region is the daughter
        params.daughters.resize(params.tracker->num_volumes());
        params.daughters.back() = std::move(daughter);

        auto u = std::make_shared<Universe>(std::move(params));
        this->insert_universe(u, std::move(built.md));
    }
}

//---------------------------------------------------------------------------//

TEST_F(GrandpaCoincidentTest, tracking)
{
    {
        // Outside heading in
        TrackResult              result = this->track({-5, 0, 0}, {1, 0, 0});
        static const char* const expected_cells[]
            = {"grandpa:outer", "baby:inner", "grandpa:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "---", "inner.s", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', ' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {4, 2, 9};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Inside heading out
        TrackResult              result = this->track({0, 0, 0}, {-1, 0, 0});
        static const char* const expected_cells[]
            = {"baby:inner", "grandpa:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "inner.s", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {1, 9};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        SCOPED_TRACE("Start on coincident point");
        TrackResult              result = this->track({-1, 0, 0}, {1, 0, 0});
        static const char* const expected_cells[]
            = {"baby:inner", "grandpa:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "inner.s", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {2, 9};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        SCOPED_TRACE("Start on coincident point (opposite direction)");
        TrackResult result = this->track({-1, 0, 0}, {-1, 0, 0});

        static const char* const expected_cells[] = {"grandpa:outer"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[] = {"---", "outer.s"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {9};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
        static const real_type expected_normals[] = {0, 0, 0, -1, 0, 0};
        EXPECT_VEC_SOFT_EQ(expected_normals, result.normals);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Simple corner case.
 *
 * This is the case for a perfect, well-connected geometry that has difficult
 * corners. It test the logic that tries to "bump" through one or more surfaces
 * at the global level while preserving the exact distance to travel by setting
 * 'movement' to a negative value.
 *
 * It's got the left corners of a flat-top hexagon, but only a single vertical
 * surface (so it's a triangular prism pointing left).
 */
class CornerTest : public TrackingGeometryTest
{
  protected:
    static constexpr real_type half_sqrt_three = 0.5 * constants::sqrt_three;
    static constexpr real_type half            = 0.5;
    static constexpr real_type apothem         = 3.0;

    void build_universes() final;
};

void CornerTest::build_universes()
{
    // Use very large tolerance to make bumps obvious
    celeritas::fuzziness() = Fuzziness{1e-2};

    ShapeContainer shapes;
    auto           ul = shapes.emplace<PlaneShape>(
        ORANGE_MD_FROM_SOURCE("ul"),
        Transform{},
        Real3{-half_sqrt_three, half, 0},
        Real3{-half_sqrt_three * apothem, half * apothem, 0});
    auto ll = shapes.emplace<PlaneShape>(
        ORANGE_MD_FROM_SOURCE("ll"),
        Transform{},
        Real3{-half_sqrt_three, -half, 0},
        Real3{-half_sqrt_three * apothem, -half * apothem, 0});
    auto ur = shapes.emplace<PlaneShape>(
        ORANGE_MD_FROM_SOURCE("ur"),
        Transform{},
        Real3{half_sqrt_three, half, 0},
        Real3{half_sqrt_three * apothem, half * apothem, 0});
    auto lr = shapes.emplace<PlaneShape>(
        ORANGE_MD_FROM_SOURCE("lr"),
        Transform{},
        Real3{half_sqrt_three, -half, 0},
        Real3{half_sqrt_three * apothem, -half * apothem, 0});

    auto hex = shapes.emplace<IntersectionShape>(
        ORANGE_MD_FROM_SOURCE("hex"),
        Transform{},
        IntersectionShape::RegionVec{
            {inside, ur}, {inside, ul}, {inside, ll}, {inside, lr}});
    auto vertical = shapes.emplace<PlaneShape>(ORANGE_MD_FROM_SOURCE("vert"),
                                               Transform{},
                                               Real3{1, 0, 0},
                                               Real3{0, 0, 0});
    // Two regions
    UnitBuilder build;
    build.exterior(
        {{inside, hex}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("exterior"));
    build.region({{inside, hex}, {inside, vertical}},
                 ZOrder::media,
                 ORANGE_MD_FROM_SOURCE("lefthex"));
    build.region({{inside, hex}, {outside, vertical}},
                 ZOrder::media,
                 ORANGE_MD_FROM_SOURCE("righthex"));

    auto built = build(std::move(ORANGE_MD_FROM_SOURCE("corner")));
    EXPECT_TRUE(built.md->is_simple());

    // Construct universe params
    Universe::Params params;
    params.id      = UniverseId{0};
    params.tracker = make_unique<SimpleUnitTracker>(std::move(built.surfaces),
                                                    std::move(built.regions));

    auto u = std::make_shared<Universe>(std::move(params));
    this->insert_universe(u, std::move(built.md));
}

//---------------------------------------------------------------------------//

TEST_F(CornerTest, tracking)
{
    // Distance to the left.corner from the center
    const real_type circum = apothem / half_sqrt_three;

    {
        // Start far left of corner, relying on bump to get us across
        auto result = this->track({-5, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"---", "corner:lefthex", "corner:righthex"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "vert.p4", "hex.p0"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[]
            = {5 - circum, circum, circum};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Start on exact left corner
        auto result = this->track({-circum, 0, 0}, {1, 0, 0});

        static const char* const expected_cells[]
            = {"corner:lefthex", "corner:righthex"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "vert.p4", "hex.p0"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '+', '+'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {circum, circum};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Start far right of corner moving left
        auto result = this->track({5, 0, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"---", "corner:righthex", "corner:lefthex"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "vert.p4", "hex.p1"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[]
            = {5 - circum, circum, circum};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }

    {
        // Start on exact right corner moving left
        auto result = this->track({circum, 0, 0}, {-1, 0, 0});

        static const char* const expected_cells[]
            = {"corner:righthex", "corner:lefthex"};
        EXPECT_VEC_EQ(expected_cells, result.cells);
        static const char* const expected_surfaces[]
            = {"---", "vert.p4", "hex.p1"};
        EXPECT_VEC_EQ(expected_surfaces, result.surfaces);
        static const char expected_senses[] = {' ', '-', '-'};
        EXPECT_VEC_EQ(expected_senses, result.senses);
        static const real_type expected_distances[] = {circum, circum};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}
