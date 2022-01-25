//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstUnitBuilder.cc
 * \brief Tests for class UnitBuilder
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../UnitBuilder.hh"

#include <iomanip>
#include "celeritas_test.hh"
#include "base/FixedViewArray.hh"
#include "base/Provenance.hh"
#include "../CuboidShape.hh"
#include "../CylinderShape.hh"
#include "../PlaneShape.hh"
#include "../SphereShape.hh"
#include "../PlacedShape.hh"
#include "../ShapeContainer.hh"

using namespace celeritas;

constexpr real_type inf = HUGE_VAL;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class UnitBuilderTest : public ::Test
{
  protected:
    //// TYPE ALIASES ////
    using SPConstShape      = UnitBuilder::SPConstShape;
    using SPConstProvenance = std::shared_ptr<const Provenance>;
    using logic_int         = CSGCell::logic_int;
    using VecSenseShape     = UnitBuilder::VecSenseShape;

    enum LogicToken : CSGCell::logic_int
    {
        LOGIC_NOT = CSGCell::LOGIC_NOT,
        LOGIC_AND = CSGCell::LOGIC_AND,
        LOGIC_OR  = CSGCell::LOGIC_OR,
    };

  protected:
    SPConstShape add_keno_cyl(const ObjectMetadata& md,
                              real_type             r,
                              real_type             zt,
                              real_type             zb)
    {
        return shapes.emplace<CylinderShape>(
            md, Transform{}, Axis::z, r, zb, zt);
    }

    SPConstShape
    add_keno_plane(const ObjectMetadata& md, const Real3& abc, real_type d)
    {
        return shapes.emplace<PlaneShape>(
            md, Transform{}, PlaneShape::from_displacement(-abc, d));
    }

  protected:
    //// DATA ////
    ShapeContainer shapes;
};

std::vector<std::string> get_surface_names(const celeritas::UnitMetadata& md)
{
    std::vector<std::string> surface_names;
    for (auto sid : range(md.num_surfaces()))
    {
        for (const auto& surface_md : md.surface_md(SurfaceId{sid}))
        {
            std::ostringstream os;
            os << sid << ":" << surface_md.first->name() << '.'
               << surface_md.second;
            surface_names.push_back(os.str());
        }
    }
    return surface_names;
}

std::ostream& operator<<(std::ostream& os, const UnitRegion& reg)
{
    using make_fixed_view;

    os << "EXPECT_VEC_EQ(CSGCell::from_string(\"" << reg.interior
       << "\").logic(), reg.interior.logic());\n"
       << "EXPECT_VEC_SOFT_EQ(Real3(" << make_fixed_view(reg.bbox.lower())
       << "), reg.bbox.lower());\n"
       << "EXPECT_VEC_SOFT_EQ(Real3(" << make_fixed_view(reg.bbox.upper())
       << "), reg.bbox.upper());\n";
    return os;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(UnitBuilderTest, functionality)
{
    // Setup
    shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("inner"), Transform(), 1.0);
    shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("outer"), Transform(), 2.0);
    shapes.emplace<CuboidShape>(ORANGE_MD_FROM_SOURCE("lcube"),
                                Transform(),
                                Real3(-.8, -.8, -.8),
                                Real3(0, .8, .8));
    shapes.emplace<CuboidShape>(ORANGE_MD_FROM_SOURCE("rcube"),
                                Transform(),
                                Real3(0, -.8, -.8),
                                Real3(.8, .8, .8));
    ASSERT_EQ(4, shapes.size());

    // Construct
    UnitBuilder build_unit;
    build_unit.exterior({{neg, shapes["outer"]}},
                        ZOrder::exterior,
                        ORANGE_MD_FROM_SOURCE("exterior"));
    build_unit.region(
        {{neg, shapes["lcube"]}}, ZOrder::hole, ORANGE_MD_FROM_SOURCE("lhole"));
    build_unit.region(
        {{neg, shapes["rcube"]}}, ZOrder::hole, ORANGE_MD_FROM_SOURCE("rhole"));
    build_unit.region({{neg, shapes["inner"]}},
                      ZOrder::media,
                      ORANGE_MD_FROM_SOURCE("interior"));
    build_unit.region({{pos, shapes["inner"]}, {neg, shapes["outer"]}},
                      ZOrder::media,
                      ORANGE_MD_FROM_SOURCE("shell"));

    auto built = build_unit(ORANGE_MD_FROM_SOURCE("unit"));

    // Test generated surfaces
    EXPECT_EQ(9, built.surfaces.size());

    // Test generated MD
    ASSERT_SP_TRUE(built.md);
    {
        const celeritas::UnitMetadata& unit_md = *built.md;
        EXPECT_EQ("unit", unit_md.metadata().name());

        // Test bounding box
        const auto& bbox = unit_md.bbox();
        EXPECT_VEC_SOFT_EQ(Real3(-2, -2, -2), bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3(2, 2, 2), bbox.upper());

        // Test cell metadata
        ASSERT_EQ(5, unit_md.num_volumes());
        {
            EXPECT_EQ("exterior", unit_md.vol_md(VolumeId(0)).name());
            EXPECT_EQ("lhole", unit_md.vol_md(VolumeId(1)).name());
            EXPECT_EQ("rhole", unit_md.vol_md(VolumeId(2)).name());
            EXPECT_EQ("interior", unit_md.vol_md(VolumeId(3)).name());
            EXPECT_EQ("shell", unit_md.vol_md(VolumeId(4)).name());
        }

        // Test surface metadata
        ASSERT_EQ(9, unit_md.num_surfaces());
        {
            const auto& surface_md = unit_md.surface_md(SurfaceId(0));
            ASSERT_EQ(1, surface_md.size());
            EXPECT_EQ(shapes["outer"], surface_md[0].first);
            EXPECT_EQ("s", surface_md[0].second);
        }
        {
            const auto& surface_md = unit_md.surface_md(SurfaceId(1));
            ASSERT_EQ(1, surface_md.size());
            EXPECT_EQ(shapes["lcube"], surface_md[0].first);
            EXPECT_EQ("mx", surface_md[0].second);
        }
        {
            const auto& surface_md = unit_md.surface_md(SurfaceId(2));
            ASSERT_EQ(2, surface_md.size());
            EXPECT_EQ(shapes["lcube"], surface_md[0].first);
            EXPECT_EQ("px", surface_md[0].second);
            EXPECT_EQ(shapes["rcube"], surface_md[1].first);
            EXPECT_EQ("mx", surface_md[1].second);
        }
    }

    // Test built regions
    ASSERT_EQ(5, built.regions.size());
    {
        const auto& reg = built.regions[0];
        EXPECT_EQ(CSGCell::from_string("0").logic(), reg.interior.logic());
        EXPECT_EQ(ZOrder::exterior, reg.zorder);
        EXPECT_EQ(geometria::infinite_bbox(), reg.bbox);
        EXPECT_FALSE(reg.has_internal_surfaces);
    }
    {
        const auto& reg = built.regions[1];
        EXPECT_EQ(CSGCell::from_string("1 2 ~ & 3 & 4 ~ & 5 & 6 ~ &").logic(),
                  reg.interior.logic());
        EXPECT_EQ(ZOrder::hole, reg.zorder);
        EXPECT_VEC_SOFT_EQ(Real3(-.8, -.8, -.8), reg.bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3(0, .8, .8), reg.bbox.upper());
        EXPECT_FALSE(reg.has_internal_surfaces);
    }
    {
        const auto& reg = built.regions[2];
        EXPECT_EQ(CSGCell::from_string("2 3 & 4 ~ & 5 & 6 ~ & 7 ~ &").logic(),
                  reg.interior.logic());
        EXPECT_EQ(ZOrder::hole, reg.zorder);
        EXPECT_VEC_SOFT_EQ(Real3(0, -.8, -.8), reg.bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3(.8, .8, .8), reg.bbox.upper());
        EXPECT_FALSE(reg.has_internal_surfaces);
    }
    {
        const auto& reg = built.regions[3];
        EXPECT_EQ(CSGCell::from_string("8 ~").logic(), reg.interior.logic());
        EXPECT_EQ(ZOrder::media, reg.zorder);
        EXPECT_VEC_SOFT_EQ(Real3(-1, -1, -1), reg.bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3(1, 1, 1), reg.bbox.upper());
        EXPECT_FALSE(reg.has_internal_surfaces);
    }
    {
        const auto& reg = built.regions[4];
        EXPECT_EQ(CSGCell::from_string("0 ~ 8 &").logic(),
                  reg.interior.logic());
        EXPECT_EQ(ZOrder::media, reg.zorder);
        EXPECT_VEC_SOFT_EQ(Real3(-2, -2, -2), reg.bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3(2, 2, 2), reg.bbox.upper());
        EXPECT_TRUE(reg.has_internal_surfaces);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Test of creating a cuboid from arbitrary planes with KENO6 input.
 *
 * This should be a reproduction of the `planes.inp` KenoVI test.
 */
TEST_F(UnitBuilderTest, keno6_planes)
{
    std::vector<SPConstShape> planes = {
        add_keno_plane(ORANGE_MD_FROM_SOURCE("mx"), {1, 0, 0}, 0.5),
        add_keno_plane(ORANGE_MD_FROM_SOURCE("px"), {-10, 0, 0}, 5.),
        add_keno_plane(ORANGE_MD_FROM_SOURCE("my"), {0, 1, 0}, 0.),
        add_keno_plane(ORANGE_MD_FROM_SOURCE("py"), {0, -1, 0}, 2.),
        add_keno_plane(ORANGE_MD_FROM_SOURCE("mz"), {0, 0, 2}, 0),
        add_keno_plane(ORANGE_MD_FROM_SOURCE("pz"), {0, 0, -1}, 3.),
    };

    VecSenseShape interior;
    for (const SPConstShape& shape : planes)
    {
        interior.emplace_back(inside, shape);
    }

    UnitBuilder build_unit;
    build_unit.exterior(
        interior, ZOrder::exterior, ORANGE_MD_FROM_SOURCE("exterior"));
    build_unit.region(
        interior, ZOrder::media, ORANGE_MD_FROM_SOURCE("interior"));

    auto built = build_unit(ORANGE_MD_FROM_SOURCE("unit"));

    const auto& bbox = built.md->bbox();
    EXPECT_VEC_SOFT_EQ(Real3(-.5, 0, 0), bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3(.5, 2, 3), bbox.upper());

    ASSERT_EQ(2, built.regions.size());
    {
        const UnitRegion& reg = built.regions[0];
        EXPECT_VEC_EQ(
            CSGCell::from_string("0 1 ~ & 2 & 3 ~ & 4 & 5 ~ & ~ ").logic(),
            reg.interior.logic())
            << "Actual: " << reg.interior;
    }
}

//---------------------------------------------------------------------------//
/*!
 * Create a single well-connected global cuboid unit.
 */
TEST_F(UnitBuilderTest, cuboid)
{
    auto cuboid = shapes.emplace<CuboidShape>(ORANGE_MD_FROM_SOURCE("cuboid"),
                                              Transform(),
                                              Real3(-.5, 0, 0),
                                              Real3(.5, 2, 3));
    UnitBuilder build_unit;
    build_unit.exterior(
        {{inside, cuboid}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("exterior"));
    build_unit.region(
        {{inside, cuboid}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("interior"));

    auto built = build_unit(ORANGE_MD_FROM_SOURCE("unit"));
    EXPECT_TRUE(built.md->is_simple());

    const auto& bbox = built.md->bbox();
    EXPECT_VEC_SOFT_EQ(Real3(-.5, 0, 0), bbox.lower());
    EXPECT_VEC_SOFT_EQ(Real3(.5, 2, 3), bbox.upper());

    ASSERT_EQ(2, built.regions.size());
    {
        // exterior should have internal surfaces
        const UnitRegion& reg = built.regions[0];
        EXPECT_VEC_EQ(
            CSGCell::from_string("0 1 ~ & 2 & 3 ~ & 4 & 5 ~ & ~ ").logic(),
            reg.interior.logic());
        EXPECT_TRUE(reg.has_internal_surfaces);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Test Matt Jessee's failing problem with explicit planes.
 *
 * The signs on the planes should be deduplicated automatically..
 */
TEST_F(UnitBuilderTest, moonraker_planes)
{
    auto c1  = add_keno_cyl(ORANGE_MD_FROM_SOURCE("1"), 4.01408e-01, 6, 0);
    auto p3  = add_keno_plane(ORANGE_MD_FROM_SOURCE("3"), {0, 0, 1}, 0.);
    auto p4  = add_keno_plane(ORANGE_MD_FROM_SOURCE("4"), {0, 0, 1}, -1.);
    auto p5  = add_keno_plane(ORANGE_MD_FROM_SOURCE("5"), {0, 0, -1}, 1.);
    auto p6  = add_keno_plane(ORANGE_MD_FROM_SOURCE("6"), {0, 0, 1}, -3.);
    auto p7  = add_keno_plane(ORANGE_MD_FROM_SOURCE("7"), {0, 0, -1}, 3.);
    auto p9  = add_keno_plane(ORANGE_MD_FROM_SOURCE("9"), {0, 0, -1}, 6.);
    auto c10 = add_keno_cyl(ORANGE_MD_FROM_SOURCE("10"), 4.098164e-01, 6, 0);
    auto c11 = add_keno_cyl(ORANGE_MD_FROM_SOURCE("11"), 4.670444e-01, 6, 0);

    UnitBuilder build_unit;
    build_unit.exterior({{inside, c11}, {inside, p3}, {inside, p9}},
                        ZOrder::exterior,
                        ORANGE_MD_FROM_SOURCE("exterior"));

    build_unit.region({{inside, c1}, {inside, p3}, {inside, p5}},
                      ZOrder::media,
                      ORANGE_MD_FROM_SOURCE("fuel.0"));
    build_unit.region({{inside, c1}, {inside, p4}, {inside, p7}},
                      ZOrder::media,
                      ORANGE_MD_FROM_SOURCE("fuel.1"));
    build_unit.region({{inside, c1}, {inside, p6}, {inside, p9}},
                      ZOrder::media,
                      ORANGE_MD_FROM_SOURCE("fuel.2"));

    build_unit.region(
        {{inside, c10}, {outside, c1}, {inside, p3}, {inside, p5}},
        ZOrder::media,
        ORANGE_MD_FROM_SOURCE("gap.3"));
    build_unit.region(
        {{inside, c10}, {outside, c1}, {inside, p4}, {inside, p7}},
        ZOrder::media,
        ORANGE_MD_FROM_SOURCE("gap.4"));
    build_unit.region(
        {{inside, c10}, {outside, c1}, {inside, p6}, {inside, p9}},
        ZOrder::media,
        ORANGE_MD_FROM_SOURCE("gap.5"));

    build_unit.region(
        {{inside, c11}, {outside, c10}, {inside, p3}, {inside, p5}},
        ZOrder::media,
        ORANGE_MD_FROM_SOURCE("clad.3"));
    build_unit.region(
        {{inside, c11}, {outside, c10}, {inside, p4}, {inside, p7}},
        ZOrder::media,
        ORANGE_MD_FROM_SOURCE("clad.4"));
    build_unit.region(
        {{inside, c11}, {outside, c10}, {inside, p6}, {inside, p9}},
        ZOrder::media,
        ORANGE_MD_FROM_SOURCE("clad.5"));

    auto built = build_unit(ORANGE_MD_FROM_SOURCE("unit1"));

    // Test surfaces
    {
        auto surface_names = get_surface_names(*built.md);
        // PRINT_EXPECTED(surface_names);
        static const char* expected_surface_names[] = {"0:11.coz",
                                                       "1:1.mz",
                                                       "1:10.mz",
                                                       "1:11.mz",
                                                       "1:3.p0",
                                                       "2:1.pz",
                                                       "2:10.pz",
                                                       "2:11.pz",
                                                       "2:9.p1",
                                                       "3:1.coz",
                                                       "4:4.p0",
                                                       "4:5.p1",
                                                       "5:6.p0",
                                                       "5:7.p1",
                                                       "6:10.coz"};
        EXPECT_VEC_EQ(expected_surface_names, surface_names);
    }

    // Test built regions
    ASSERT_EQ(10, built.regions.size());
    {
        const UnitRegion& reg = built.regions[0];
        EXPECT_VEC_EQ(CSGCell::from_string("1 0 ~ 1 & 2 ~ & & 2 ~ & ~").logic(),
                      reg.interior.logic())
            << "Actual: " << reg.interior;
        EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), reg.bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), reg.bbox.upper());
    }
}

//---------------------------------------------------------------------------//
/*!
 * Test *shape* deduplication and implicit exterior.
 */
TEST_F(UnitBuilderTest, duplicate)
{
    auto make_cube = [this](ObjectMetadata md) {
        return this->shapes.emplace<CuboidShape>(
            std::move(md), Transform(), Real3(-1, -1, -1), Real3(1, 1, 1));
    };

    auto intcube = make_cube(ORANGE_MD_FROM_SOURCE("interior"));
    auto extcube = make_cube(ORANGE_MD_FROM_SOURCE("exterior"));

    UnitBuilder build_unit;
    build_unit.exterior({{inside, extcube}},
                        ZOrder::implicit_exterior,
                        ORANGE_MD_FROM_SOURCE("boundary"));
    build_unit.region(
        {{inside, intcube}}, ZOrder::media, ORANGE_MD_FROM_SOURCE("media"));

    auto built = build_unit(ORANGE_MD_FROM_SOURCE("unit"));

    {
        auto surface_names = get_surface_names(*built.md);
        // PRINT_EXPECTED(surface_names);
        static const char* expected_surface_names[] = {"0:exterior.mx",
                                                       "0:interior.mx",
                                                       "1:exterior.px",
                                                       "1:interior.px",
                                                       "2:exterior.my",
                                                       "2:interior.my",
                                                       "3:exterior.py",
                                                       "3:interior.py",
                                                       "4:exterior.mz",
                                                       "4:interior.mz",
                                                       "5:exterior.pz",
                                                       "5:interior.pz"};
        EXPECT_VEC_EQ(expected_surface_names, surface_names);
        EXPECT_TRUE(built.md->is_simple());
    }
    {
        const UnitRegion& reg = built.regions[0];
        EXPECT_VEC_EQ(CSGCell::from_string("* ~").logic(), reg.interior.logic())
            << "Actual: " << reg.interior;
        EXPECT_VEC_SOFT_EQ(Real3({-inf, -inf, -inf}), reg.bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({inf, inf, inf}), reg.bbox.upper());
    }
    {
        const UnitRegion& reg = built.regions[1];
        EXPECT_VEC_EQ(CSGCell::from_string("*").logic(), reg.interior.logic())
            << "Actual: " << reg.interior;
        EXPECT_VEC_SOFT_EQ(Real3({-1, -1, -1}), reg.bbox.lower());
        EXPECT_VEC_SOFT_EQ(Real3({1, 1, 1}), reg.bbox.upper());
        EXPECT_FALSE(reg.has_internal_surfaces);
    }
}
