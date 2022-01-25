//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/XMLUtils.hh
 * \brief XMLUtils class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include <string>
#include "base/Provenance.hh"
#include "orange/Definitions.hh"
#include "orange/query/ObjectMetadata.hh"

namespace Teuchos
{
class ParameterList;
}

namespace celeritas
{
class PlacedShape;
namespace detail
{
//---------------------------------------------------------------------------//
// Build metadata from a plist (name, optional description, provenance)
ObjectMetadata build_md(const Teuchos::ParameterList& plist);

// Build a shape from the given plist
std::shared_ptr<const PlacedShape> build_shape(const Teuchos::ParameterList&);

// Get a 3-vector of doubles from a plist item
Real3 get_space_vector(const Teuchos::ParameterList&, const std::string& key);

// Get a 3-vector of indices from a plist item
Array<def::size_type, 3>
get_dim_vector(const Teuchos::ParameterList& plist, const std::string& key);

// Get a 3x3 matrix from a plist item
SpaceMatrix
get_space_matrix(const Teuchos::ParameterList&, const std::string& key);

// Construct a transform from translate/rotate
Transform build_transform(const Teuchos::ParameterList&);

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
