//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/test/tstPolygonUtils.cc
 * \brief PolygonUtils unit tests
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../PolygonUtils.hh"

#include "celeritas_test.hh"
#include "base/Constants.hh"
#include "orange/surfaces/SurfaceContainer.hh"
#include "orange/Fuzziness.hh"

using celeritas::detail::edges_intersect;
using celeritas::detail::edges_overlap;
using celeritas::detail::is_right_turn;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class PolygonUtilsTest : public Test
{
  protected:
    //// TYPEDEFS ////
    using Point    = Array<real_type, 2>;
    using VecPoint = std::vector<Point>;

    void SetUp() override
    {
        celeritas::fuzziness() = celeritas::Fuzziness{1e-6};
    }

    VecPoint create_rhexagon(real_type apothem, const Point& center) const;
};

//---------------------------------------------------------------------------//
/*!
 * Create a polygon representation of a rotated hexagon
 *
 * \verbatim
 *   _____
 *  /     \
 * /       \x   The vertices are ordered counter clockwise
 * \       /    starting at the +x axis vertex denoted x.
 *  \_____/
 * \endverbatim
 */

auto PolygonUtilsTest::create_rhexagon(real_type    apothem,
                                       const Point& center) const -> VecPoint
{
    CELER_EXPECT(apothem >= 0);
    using constants::pi;
    int       num_sides = 6;
    real_type r         = apothem / std::cos(pi / num_sides);

    // Create polygon points.
    return {center + Point{r, 0},
            center + Point{half * r, apothem},
            center + Point{-half * r, apothem},
            center + Point{-r, 0},
            center + Point{-half * r, -apothem},
            center + Point{half * r, -apothem}};
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(PolygonUtilsTest, edges_overlap)
{
    Point point_A(1, 0);
    Point point_B(2, 0);
    Point point_C(3, 0);

    EXPECT_FALSE(edges_overlap(point_A, point_B, point_C));
    EXPECT_TRUE(edges_overlap(point_A, point_B, point_A));
}

TEST_F(PolygonUtilsTest, edges_intersect)
{
    // Same line test
    {
        Point point_A(1, 0);
        Point point_B(2, 0);
        Point point_C(2, 0);
        Point point_D(2, 1);
        EXPECT_TRUE(edges_intersect(point_A, point_B, point_A, point_B));
        EXPECT_TRUE(edges_intersect(point_A, point_B, point_B, point_A));
    }

    // Parallel non-intersecting line test
    {
        Point point_A(1, 0);
        Point point_B(2, 0);
        Point point_C(1, 1);
        Point point_D(2, 1);
        EXPECT_FALSE(edges_intersect(point_A, point_B, point_C, point_D));
    }

    // Lines do Not intersect test
    {
        Point point_A(1, 0);
        Point point_B(2, 0);
        Point point_C(2, 1);
        Point point_D(2, 2);
        EXPECT_FALSE(edges_intersect(point_A, point_B, point_C, point_D));
        EXPECT_FALSE(edges_intersect(point_C, point_D, point_A, point_B));
    }

    // Lines do intersect test
    {
        Point point_A(1, 1);
        Point point_B(3, 1);
        Point point_C(2, 0);
        Point point_D(2, 2);
        Point point_E(2, 1);
        Point point_F(2, 2);
        EXPECT_TRUE(edges_intersect(point_A, point_B, point_C, point_D));
        EXPECT_TRUE(edges_intersect(point_A, point_B, point_E, point_F));
    }
}

TEST_F(PolygonUtilsTest, is_right_turn)
{
    Point point_A(1, 0);
    Point point_B1(2, 0);
    Point point_B2(2, -1);
    Point point_B3(2, 1);
    Point point_C(3, 0);

    EXPECT_FALSE(is_right_turn(point_A, point_B1, point_C));
    EXPECT_FALSE(is_right_turn(point_A, point_B2, point_C));
    EXPECT_TRUE(is_right_turn(point_A, point_B3, point_C));
}
