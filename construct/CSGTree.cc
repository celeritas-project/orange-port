//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CSGTree.cc
 * \brief CSGTree class definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CSGTree.hh"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include "base/ContainerFunctions.hh"
#include "base/Range.hh"
#include "CSGCell.hh"

namespace celeritas
{
namespace
{
struct LessDaughter
{
    using Daughter = CSGTree::Daughter;

    bool operator()(const Daughter& lhs, const Daughter& rhs) const
    {
        return std::make_tuple(lhs.second, lhs.first)
               < std::make_tuple(rhs.second, rhs.first);
    }
};
} // namespace

//---------------------------------------------------------------------------//
/*!
 * Insert a node with any number of daughters.
 *
 * An empty node is always true (infinite) no matter the conjunction.
 *
 * TODO: we eventually normalize (!A | !B) -> !(A & B) for node deduplication?
 * Probably after adding a "not" node.
 */
auto CSGTree::emplace(VecDaughter daughters, LogicToken conj) -> NodeId
{
    CELER_EXPECT(std::all_of(
        daughters.begin(), daughters.end(), [this](const Daughter& d) {
            return d.second < this->size();
        }));

    return this->insert_impl(this->build_node(std::move(daughters), conj));
}

//---------------------------------------------------------------------------//
/*!
 * Replace a node (and its daughters if applicable) with a logical value.
 */
void CSGTree::replace(NodeId node, bool replacement)
{
    CELER_EXPECT(node < this->size());

    // Recursively replace daughters as needed
    size_type min_id = this->replace_impl(node.get(), replacement);

    // Simplify all IDs higher than the minimum
    for (auto i : range(min_id, this->size()))
    {
        const CSGNode& n = nodes_[i];
        if (!n.is_leaf())
        {
            nodes_[i] = this->build_node(n.daughters(), n.conjunction());
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * Convert a node to a CSG cell expression.
 */
CSGCell CSGTree::build_cell(NodeId i) const
{
    CELER_EXPECT(i < this->size());
    VecLogic logic;
    this->build_node_logic(i, &logic);
    CELER_ENSURE(!logic.empty());
    return CSGCell(std::move(logic));
}

//---------------------------------------------------------------------------//
/*!
 * Describe all nodes.
 */
void CSGTree::describe(std::ostream& os) const
{
    using std::left;
    using std::setw;

    std::string sep(4 + 1 + 60 + 1, '=');
    sep[4]     = ' ';
    sep.back() = '\n';

    os << sep;
    for (auto i : range(this->size()))
    {
        os << setw(4) << i << ' ';

        const CSGNode& node = this->at(NodeId{i});
        if (node.is_surface())
        {
            os << "Surface " << node.surface_id().get();
        }
        else if (node.is_true())
        {
            os << "[true]";
        }
        else
        {
            CELER_ASSERT(!node.is_leaf());
            bool first_entry = true;
            char conj        = (node.conjunction() == CSGCell::LOGIC_AND  ? '&'
                                : node.conjunction() == CSGCell::LOGIC_OR ? '|'
                                                                          : '?');

            for (const auto& daughter : node.daughters())
            {
                if (!first_entry)
                {
                    os << conj;
                }
                else
                {
                    first_entry = false;
                }

                if (daughter.first == inside)
                {
                    os << '~';
                }
                os << daughter.second.get();
            }
        }
        os << '\n';
    }
    os << sep << '\n';
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
/*!
 * Create a node with multiple daughters.
 */
auto CSGTree::build_node(VecDaughter daughters, LogicToken conj) -> CSGNode
{
    if (daughters.empty())
    {
        // Everywhere; exit early
        return CSGNode{};
    }

    // Move any constants to the end of the vector for processing
    auto cdaughters = std::partition(
        daughters.begin(), daughters.end(), [this](const Daughter& d) {
            return !this->at(d.second).is_true();
        });

    // Check for constants that simplify
    for (auto iter = cdaughters; iter != daughters.end(); ++iter)
    {
        CELER_ASSERT(this->at(iter->second).is_true());

        if (conj == CSGCell::LOGIC_OR && iter->first == Sense::inside)
        {
            // Inside any or "everywhere" -> everywhere; exit early
            return CSGNode{};
        }
        else if (conj == CSGCell::LOGIC_AND && iter->first == Sense::outside)
        {
            // Inside any and nowhere -> all false; exit early
            auto everywhere = this->insert_impl(CSGNode{});
            return CSGNode({outside, everywhere});
        }
    }

    if (cdaughters == daughters.begin())
    {
        // Only constants are left, but there's no "shortcut". That means all
        // the constants are the same and singular.
        return this->build_node(*cdaughters);
    }

    // Only useless constants remain: delete them
    daughters.erase(cdaughters, daughters.end());
    CELER_ASSERT(!daughters.empty());

    // Sort remaining surfaces/nodes by increasing ID and eliminate duplicates
    std::sort(daughters.begin(), daughters.end(), LessDaughter{});
    uniquify(daughters, LessDaughter{});

    if (daughters.size() == 1)
    {
        // Jump to the version that will simplify single-element daughters
        return this->build_node(daughters.front());
    }

    return CSGNode{daughters, conj};
}

//---------------------------------------------------------------------------//
/*!
 * Insert a single daughter.
 */
CSGNode CSGTree::build_node(Daughter d)
{
    CELER_EXPECT(d.second < this->size());

    // See if daughter is a node that points to a leaf; if so, combine the
    // senses (xor) and point to the leaf ourselves.
    const CSGNode& dnode = this->at(d.second);
    if (d.first == inside)
    {
        // Identical to daughter
        return dnode;
    }
    else if (!dnode.is_leaf() && dnode.daughters().size() == 1)
    {
        // Single-daughter node can only be inverse, so we're the same as the
        // granddaughter node
        // TODO: simplify silly "not" logic
        const Daughter& granddaughter = dnode.daughters().front();
        CELER_ASSERT(granddaughter.first == outside);
        return this->at(granddaughter.second);
    }

    return CSGNode{d};
}

//---------------------------------------------------------------------------//
/*!
 * Insert a node, with deduplication.
 */
auto CSGTree::insert_impl(CSGNode n) -> NodeId
{
    // Insert and search for existing node
    auto iter_inserted = ids_.emplace(std::move(n), NodeId(this->size()));

    if (iter_inserted.second)
    {
        // Copy newly created CSG node to list
        nodes_.push_back(iter_inserted.first->first);
    }

    CELER_ENSURE(iter_inserted.first->second < this->size());
    return iter_inserted.first->second;
}

//---------------------------------------------------------------------------//
/*!
 * Recursively propagate "known true" or "known false" nodes to daughters.
 *
 * Return the lowest node ID encountered
 */
auto CSGTree::replace_impl(size_type i, bool replacement) -> size_type
{
    size_type min_id = i;
    CSGNode   orig   = nodes_[i];

    LogicToken replace_if_conj;
    Sense      replace_sense;
    if (replacement)
    {
        nodes_[i] = CSGNode{};

        // 'always true' with 'and' implies that all "inside" daughters are
        // true and "outside" daughters are false.
        replace_if_conj = CSGCell::LOGIC_AND;
        replace_sense   = inside;
    }
    else
    {
        auto everywhere = this->insert_impl(CSGNode{});
        nodes_[i]       = CSGNode({outside, everywhere});

        // 'always true' with 'and' implies that all "inside" daughters are
        // true and "outside" daughters are false.
        replace_if_conj = CSGCell::LOGIC_OR;
        replace_sense   = outside;
    }

    if (!orig.is_leaf()
        && (orig.conjunction() == replace_if_conj
            || orig.daughters().size() == 1))
    {
        // Replace daughters based on logical induction
        for (const auto& d : orig.daughters())
        {
            min_id = std::min(
                min_id,
                this->replace_impl(d.second.get(), d.first == replace_sense));
        }
    }

    return min_id;
}

//---------------------------------------------------------------------------//
/*!
 * Recursively construct RPN logic for a node.
 */
void CSGTree::build_node_logic(NodeId i, VecLogic* logic) const
{
    CELER_EXPECT(logic);
    const CSGNode& node = this->at(i);
    if (node.is_surface())
    {
        // This "node" is just a single surface
        logic->push_back(node.surface_id().get());
    }
    else if (node.is_true())
    {
        // This "node" is a constant. Should not be called if the tree is fully
        // simplified, except perhaps at the global level.
        logic->push_back(CSGCell::LOGIC_TRUE);
    }
    else
    {
        CELER_ASSERT(!node.is_leaf());
        // Process daughters
        bool       first_entry = true;
        LogicToken conjunction = node.conjunction();

        for (const auto& daughter : node.daughters())
        {
            // Recurse into daughter
            this->build_node_logic(daughter.second, logic);

            // Add sense for this daughter
            if (daughter.first == outside)
            {
                if (logic->empty() || logic->back() != CSGCell::LOGIC_NOT)
                {
                    logic->push_back(CSGCell::LOGIC_NOT);
                }
                else
                {
                    // Reverse trailing not, possibly from a single surface
                    logic->pop_back();
                }
            }

            if (!first_entry)
            {
                logic->push_back(conjunction);
            }

            first_entry = false;
        }
    }
}
//---------------------------------------------------------------------------//
} // namespace celeritas
