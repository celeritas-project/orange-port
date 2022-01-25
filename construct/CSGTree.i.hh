//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CSGTree.i.hh
 * \brief CSGTree inline method definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Insert a single daughter.
 */
auto CSGTree::emplace(Daughter d) -> NodeId
{
    CELER_EXPECT(d.second < this->size());
    return this->insert_impl(this->build_node(d));
}

//---------------------------------------------------------------------------//
/*!
 * Access a node.
 */
const CSGNode& CSGTree::at(NodeId i) const
{
    CELER_EXPECT(i < nodes_.size());
    return nodes_[i.get()];
}

//---------------------------------------------------------------------------//
} // namespace celeritas
