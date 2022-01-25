//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/test/tstFoundCells.cc
 * \brief Tests for class FoundCells
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../FoundCells.hh"

#include "celeritas_test.hh"
#include "orange/Definitions.hh"

using celeritas::detail::FoundCells;
using namespace celeritas;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(FoundCellsTest, empty)
{
    FoundCells fv;
    EXPECT_TRUE(fv.empty());
    EXPECT_EQ(0, fv.size());
    EXPECT_EQ(zorder_int(ZOrder::invalid), fv.zorder());
    EXPECT_EQ(0, fv.cell_view().size());
}

//---------------------------------------------------------------------------//
TEST(FoundCellsTest, multi_zorder)
{
    FoundCells fv;

    // First insertion succeeds (not flipped face)
    EXPECT_TRUE(fv.insert(VolumeId{2},
                          zorder_int(ZOrder::hole),
                          FoundFace{true, neg, FaceId{3}}));
    EXPECT_EQ(1, fv.size());
    EXPECT_FALSE(fv.empty());

    // Second insertion is a lower z order and should cause cell search
    // to end
    EXPECT_FALSE(
        fv.insert(VolumeId{3}, zorder_int(ZOrder::media), FoundFace{true}));

    // Check results
    EXPECT_EQ(1, fv.size());
    EXPECT_EQ(zorder_int(ZOrder::hole), fv.zorder());
    ASSERT_EQ(1, fv.cell_view().size());
    EXPECT_EQ(2, fv.cell_view().front().get());
    EXPECT_EQ(3, fv.face().unchecked_get());
    EXPECT_EQ(Sense::neg, fv.sense());
}

//---------------------------------------------------------------------------//
TEST(FoundCellsTest, overlapping_same)
{
    FoundCells fv;

    // First insertion succeeds
    EXPECT_TRUE(fv.insert(VolumeId{1},
                          zorder_int(ZOrder::hole),
                          FoundFace{true, pos, FaceId{2}}));
    EXPECT_EQ(1, fv.size());
    EXPECT_FALSE(fv.empty());
    EXPECT_EQ(Sense::pos, fv.sense());

    // Overlapping hole (opposite side of surface)
    EXPECT_TRUE(fv.insert(VolumeId{3},
                          zorder_int(ZOrder::hole),
                          FoundFace{true, pos, FaceId{4}}));
    EXPECT_EQ(2, fv.size());
    EXPECT_FALSE(fv.empty());

    // XXX implementation detail: max size = 2, so insertion fails
    EXPECT_FALSE(fv.insert(VolumeId{5},
                           zorder_int(ZOrder::hole),
                           FoundFace{true, pos, FaceId{4}}));
    EXPECT_EQ(2, fv.size());

    // Test overlapping values
    EXPECT_EQ(1, fv.cell_view()[0].get());
    EXPECT_EQ(3, fv.cell_view()[1].get());
    EXPECT_EQ(zorder_int(ZOrder::hole), fv.zorder());
}

//---------------------------------------------------------------------------//

TEST(FoundCellsTest, errors)
{
#ifndef REQUIRE_ON
    SKIP_TEST("DBC checking is disabled");
#endif
    FoundCells fv;
    // Invalid cell ID
    EXPECT_THROW(fv.insert(VolumeId{},
                           zorder_int(ZOrder::hole),
                           FoundFace{true, pos, FaceId{1}}),
                 assertion);
    ASSERT_TRUE(fv.empty());

    // Not-found face
    EXPECT_THROW(fv.insert(VolumeId{3}, zorder_int(ZOrder::hole), FoundFace{}),
                 assertion);
    ASSERT_TRUE(fv.empty());

    fv.insert(VolumeId{4},
              zorder_int(ZOrder::media),
              FoundFace{true, pos, FaceId{4}});
    // Increasing zorder
    EXPECT_THROW(fv.insert(VolumeId{2},
                           zorder_int(ZOrder::hole),
                           FoundFace{true, pos, FaceId{2}}),
                 assertion);
}
