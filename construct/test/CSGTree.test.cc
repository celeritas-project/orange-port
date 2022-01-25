//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/tstCSGTree.cc
 * \brief Tests for class CSGTree
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../CSGTree.hh"

#include "celeritas_test.hh"

using celeritas::CSGCell;
using celeritas::CSGTree;
using celeritas::inside;
using celeritas::outside;
using celeritas::SurfaceId;

static const auto AND = celeritas::CSGCell::LOGIC_AND;
static const auto OR  = celeritas::CSGCell::LOGIC_OR;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class CSGTreeTest : public ::Test
{
  protected:
    //// TYPE ALIASES ////
    using VecInt   = std::vector<int>;
    using NodeId   = CSGTree::NodeId;
    using Daughter = CSGTree::Daughter;
    using VecLogic = CSGCell::VecLogic;

  protected:
    void SetUp() {}

    VecInt add_surfaces(VecInt ids)
    {
        VecInt result;
        // Return node ID values
        for (int i : ids)
        {
            auto node_id = tree.emplace(SurfaceId(i));
            result.push_back(node_id.get());
        }
        return result;
    }

    VecLogic to_logic(const char* s)
    {
        return CSGCell::from_string(s).logic();
    }

    VecLogic to_logic(NodeId i) { return tree.build_cell(i).logic(); }

  protected:
    //// DATA ////
    CSGTree tree;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CSGTreeTest, surfaces)
{
    EXPECT_VEC_EQ(VecInt({0, 1, 2, 3, 4}), this->add_surfaces({0, 1, 3, 4, 5}));
    EXPECT_EQ(5, tree.size());
    const auto& n = tree.at(NodeId(2));
    EXPECT_TRUE(n.is_leaf());
    EXPECT_TRUE(n.is_surface());
    EXPECT_EQ(3, n.surface_id().get());

    // Adding surfaces again should get same result
    EXPECT_VEC_EQ(VecInt({0, 1, 2, 3, 4}), this->add_surfaces({0, 1, 3, 4, 5}));
    EXPECT_EQ(5, tree.size());
}

TEST_F(CSGTreeTest, simple)
{
    this->add_surfaces({0, 1, 2, 3, 4});
    // Add single-surface daughters
    auto i1 = tree.emplace({{inside, NodeId{0}}}, AND);
    EXPECT_EQ(0, i1.get());
    auto i1a = tree.emplace(Daughter{inside, NodeId{0}});
    EXPECT_EQ(i1.get(), i1a.get());

    // Add a node that should simplify to a single surface
    auto i1b = tree.emplace(Daughter{inside, i1});
    EXPECT_EQ(i1.get(), i1b.get());

    // Add equivalent multi-node daughters
    auto i2 = tree.emplace({{inside, NodeId{0}}, {outside, NodeId{1}}}, AND);
    EXPECT_EQ(5, i2.get());
    auto i2a = tree.emplace({{outside, NodeId{1}}, {inside, NodeId{0}}}, AND);
    EXPECT_EQ(i2.get(), i2a.get());

    // Add another node
    auto i3 = tree.emplace({{inside, NodeId{2}}, {inside, NodeId{3}}}, AND);
    EXPECT_EQ(6, i3.get());

    // Add a node of nodes
    auto i4 = tree.emplace({{inside, i2}, {inside, i3}}, AND);
    EXPECT_EQ(this->to_logic("0 1 ~ & 2 3 & &"), this->to_logic(i4));

    // Add a node that's identical to the daughter
    auto i4a = tree.emplace({inside, i4});
    EXPECT_EQ(i4.get(), i4a.get());

    // *not* identical to daughter
    auto not_i4 = tree.emplace({outside, i4});
    EXPECT_EQ(i4.get() + 1, not_i4.get());

    // Double inverse: same as i4
    auto i4b = tree.emplace({outside, not_i4});
    EXPECT_EQ(i4.get(), i4b.get());

    tree.describe(cout);
}

TEST_F(CSGTreeTest, replace)
{
    this->add_surfaces({0, 1, 2, 3});
    auto o0  = tree.emplace(Daughter{outside, NodeId{0}});
    auto o1  = tree.emplace(Daughter{outside, NodeId{1}});
    auto i01 = tree.emplace({{inside, NodeId{0}}, {inside, NodeId{1}}}, AND);
    auto o01 = tree.emplace({{outside, NodeId{0}}, {outside, NodeId{1}}}, AND);
    auto o02 = tree.emplace({{outside, NodeId{0}}, {outside, NodeId{2}}}, AND);
    auto i2or = tree.emplace({{outside, NodeId{0}}, {inside, NodeId{2}}}, OR);
    auto n3   = tree.emplace({{inside, NodeId{0}}, {outside, NodeId{3}}}, AND);
    EXPECT_EQ(4 + 2 + 5, tree.size());

    // tree.describe(cout);

    // Assume *always* inside both 0 and 1
    tree.replace(i01, true);
    EXPECT_EQ(this->to_logic("*"), this->to_logic(NodeId{0}));
    EXPECT_EQ(this->to_logic("* ~"), this->to_logic(o0));
    EXPECT_EQ(this->to_logic("*"), this->to_logic(NodeId{1}));
    EXPECT_EQ(this->to_logic("* ~"), this->to_logic(o1));
    EXPECT_EQ(this->to_logic("*"), this->to_logic(i01));
    EXPECT_EQ(this->to_logic("* ~"), this->to_logic(o01));
    EXPECT_EQ(this->to_logic("* ~"), this->to_logic(o02));
    EXPECT_EQ(this->to_logic("2"), this->to_logic(i2or));
    EXPECT_EQ(this->to_logic("3 ~"), this->to_logic(n3));

    // Replace "not three" with false -> 3 becomes true
    tree.replace(n3, false);
    EXPECT_EQ(this->to_logic("*"), this->to_logic(NodeId{3}));

    // tree.describe(cout);
}

// Surrogate for UnitProtoTest.build_array_parent
TEST_F(CSGTreeTest, replace_complex) {}
