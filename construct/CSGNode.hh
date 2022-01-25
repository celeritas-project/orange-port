//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CSGNode.hh
 * \brief CSGNode class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <functional>
#include <iosfwd>
#include <utility>
#include <vector>

#include "base/OpaqueId.hh"
#include "orange/Definitions.hh"
#include "CSGCell.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Low-level data structure for constructing a CSG cell.
 *
 * A CSG node in ORANGE is designed to be friendly to the "shapes are many
 * faces" and "regions are the intersection of many shapes"
 * - A leaf, which can be either a surface or a special constant; or
 * - A node (set of node+sense daughters) joined by 'and' or 'or'.
 */
class CSGNode
{
  public:
    //@{
    //! Types
    using CSGLogicToken = CSGCell::LogicToken;
    using NodeId        = OpaqueId<CSGNode>;
    using Daughter      = std::pair<Sense, NodeId>;
    using VecDaughter   = std::vector<Daughter>;
    //@}

  public:
    // Construct a "true" leaf
    inline CSGNode() = default;

    // Construct a leaf (single surface)
    explicit inline CSGNode(SurfaceId id);

    // Construct a node (single daughters)
    explicit inline CSGNode(Daughter daughter);

    // Construct a node (multiple daughters)
    CSGNode(VecDaughter daughter, CSGCell::LogicToken conjunction);

    //// ACCESSORS ////

    //! Whether this entry is a leaf of the CSG tree
    bool is_leaf() const { return is_surface() || is_true(); }

    //! Whether this entry corresponds to a single surface
    bool is_surface() const { return type_ == NodeType::leaf_surface; }

    //! Whether this entry is the constant "true"
    bool is_true() const { return type_ == NodeType::leaf_true; }

    // Surface ID
    inline SurfaceId surface_id() const;

    // Daughter components of this node
    inline const VecDaughter& daughters() const;

    //! The conjunction type if this is a node
    inline CSGLogicToken conjunction() const;

    // Equality with another CSG node
    inline bool operator==(const CSGNode& other) const;

  private:
    using surfid_int = SurfaceId::size_type;

    enum class NodeType
    {
        leaf_true,    // Always true
        leaf_surface, // "ID" is a surface ID
        node_and,     // "ID" is the index in cells_, with 'and' conjunction
        node_or,      // Same but with 'or' conjunction
        // TODO: node_not
    };

    NodeType    type_ = NodeType::leaf_true;
    surfid_int  id_   = static_cast<surfid_int>(-1);
    VecDaughter daughters_;

    friend struct std::hash<CSGNode>;
};

// Print to a stream.
std::ostream& operator<<(std::ostream& os, const CSGNode& n);

//---------------------------------------------------------------------------//
} // namespace celeritas

// Hash specialization
namespace std
{
template<>
struct hash<celeritas::CSGNode>
{
    using argument_type = celeritas::CSGNode;
    using result_type   = size_type;

    result_type operator()(const argument_type&) const;
};
} // namespace std

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "CSGNode.i.hh"
//---------------------------------------------------------------------------//
