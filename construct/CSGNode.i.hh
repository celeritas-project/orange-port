//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CSGNode.i.hh
 * \brief CSGNode inline method definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct a leaf (single surface)
 */
CSGNode::CSGNode(SurfaceId id) : type_(NodeType::leaf_surface), id_(id.get())
{
    CELER_ENSURE(this->is_surface());
}

//---------------------------------------------------------------------------//
/*!
 * Construct a node (single surface)
 */
CSGNode::CSGNode(Daughter daughter) : CSGNode({daughter}, CSGCell::LOGIC_AND)
{
    CELER_ENSURE(!this->is_leaf());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Surface ID
 */
auto CSGNode::surface_id() const -> SurfaceId
{
    CELER_EXPECT(this->is_leaf());
    return SurfaceId{id_};
}

//---------------------------------------------------------------------------//
/*!
 * \brief The daughter elements if this is a node
 */
auto CSGNode::daughters() const -> const VecDaughter&
{
    CELER_EXPECT(!this->is_leaf());
    return daughters_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief The conjunction type if this is a node
 */
auto CSGNode::conjunction() const -> CSGLogicToken
{
    CELER_EXPECT(!this->is_leaf());
    return (type_ == NodeType::node_and ? CSGCell::LOGIC_AND
                                        : CSGCell::LOGIC_OR);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Equality with another CSG node.
 */
bool CSGNode::operator==(const CSGNode& other) const
{
    return type_ == other.type_ && id_ == other.id_
           && daughters_ == other.daughters_;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
