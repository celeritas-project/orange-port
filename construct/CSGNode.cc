//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CSGNode.cc
 * \brief CSGNode class definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CSGNode.hh"

#include <iostream>
#include "Nemesis/serialize/FNVHasher.hh"
#include "Nemesis/serialize/std_pair.hh"
#include "Nemesis/serialize/std_vector.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct a node.
 *
 * This is only used by CSGTree.
 */
CSGNode::CSGNode(VecDaughter daughters, CSGCell::LogicToken conjunction)
    : type_(conjunction == CSGCell::LOGIC_AND ? NodeType::node_and
                                              : NodeType::node_or)
    , daughters_(std::move(daughters))
{
    CELER_EXPECT(!daughters_.empty());
    CELER_EXPECT(conjunction == CSGCell::LOGIC_AND
                 || conjunction == CSGCell::LOGIC_OR);
    CELER_ENSURE(!this->is_leaf());
}

//---------------------------------------------------------------------------//
/*!
 * Print to a stream.
 */
std::ostream& operator<<(std::ostream& os, const CSGNode& n)
{
    os << "Node{";
    if (n.is_surface())
    {
        os << n.surface_id().get();
    }
    else if (!n.is_leaf())
    {
        for (const auto& d : n.daughters())
        {
            os << '{' << to_char(d.first) << d.second.get() << "},";
        }
    }
    os << '}';
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas

namespace std
{
//---------------------------------------------------------------------------//
/*!
 * Hash specialization.
 */
auto hash<celeritas::CSGNode>::operator()(const argument_type& n) const
    -> result_type
{
    FNVHasher<sizeof(result_type)> hasher;

    hasher << n.type_;
    if (n.is_leaf())
    {
        hasher << n.id_;
    }
    else
    {
        hasher << n.daughters_;
    }
    return hasher.hash();
}

} // namespace std
