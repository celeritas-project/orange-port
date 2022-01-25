//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/Transform.hh
 * \brief Transform class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "../Definitions.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
void transform_to_daughter(SpanConstReal3   parent_pos,
                           SpanConstReal3   parent_dir,
                           const Transform& transform,
                           SpanReal3        pos,
                           SpanReal3        dir);

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
