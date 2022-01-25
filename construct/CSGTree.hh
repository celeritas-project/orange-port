//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CSGTree.hh
 * \brief CSGTree class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <unordered_map>
#include <vector>
#include "CSGNode.hh"
#include "CSGCell.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Represention of a tree of constructive solid geometry elements.
 *
 * A typical "node" is a single shape, but an "IntersectionShape", created by
 * inserting holes, defining boundaries, and adding chords, is a shape of
 * shapes: a node of depth 2.
 *
 * After everything is inserted, simplifications are performed:
 * - External boundary shape is replaced with 'true' (always inside)
 * - Implications of 'true' shape are propagated to its components (e.g.
 *   turning +Z into 'false')
 * - Tree is simplified from bottom up, removing 'true' or 'false' nodes
 *   where applicable
 *
 * \todo:
 *  - add object metadata
 *  - instead of daughter being sense + node, have an other 'node type' that's
 *    simply negation, which will reduce storage and simplify construction
 *  - decompose into 'builder' (with hash table) and 'built' (transfer to unit
 *    metadata for nice output)
 */
class CSGTree
{
  public:
    //@{
    //! Public type aliases
    using NodeId      = OpaqueId<CSGNode>;
    using Daughter    = CSGNode::Daughter;
    using VecDaughter = CSGNode::VecDaughter;
    using LogicToken  = CSGCell::LogicToken;
    using size_type   = NodeId::size_type;
    //@}

  public:
    // Constructor
    CSGTree() = default;

    //! Insert an "everywhere" node
    NodeId emplace() { return this->insert_impl(CSGNode{}); }

    //! Insert a surface node
    NodeId emplace(SurfaceId id) { return this->insert_impl(CSGNode{id}); }

    // Insert a single daughter
    inline NodeId emplace(Daughter d);

    // Insert multiple daughters
    NodeId emplace(VecDaughter daughter, LogicToken conjunction);

    // Access a node
    inline const CSGNode& at(NodeId i) const;

    //! Number of nodes
    size_type size() const { return nodes_.size(); }

    //// ADVANCED ////

    // Convert a node to a CSG cell expression
    CSGCell build_cell(NodeId i) const;

    // Describe all nodes
    void describe(std::ostream& os) const;

    // Replace a node (and its daughters if applicable) with a logical value
    void replace(NodeId i, bool replacement);

  private:
    using VecLogic = CSGCell::VecLogic;

    //// DATA ////

    // Logical expressions for each completed "node" of the CSG tree
    std::vector<CSGNode> nodes_;

    // Hashed nodes (originals + simplified) -> node IDs
    std::unordered_map<CSGNode, NodeId> ids_;

    //// FUNCTIONS ////

    CSGNode build_node(Daughter d);
    CSGNode build_node(VecDaughter daughter, LogicToken conjunction);

    NodeId    insert_impl(CSGNode node);
    size_type replace_impl(size_type i, bool replacement);

    void build_node_logic(NodeId i, VecLogic* logic) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "CSGTree.i.hh"
//---------------------------------------------------------------------------//
