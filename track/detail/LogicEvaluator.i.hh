//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/LogicEvaluator.i.hh
 * \brief LogicEvaluator inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Macros.hh"
#include "base/Assert.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Construct with a view to some logic definition
 */
CELER_FORCEINLINE_FUNCTION LogicEvaluator::LogicEvaluator(SpanConstLogic logic)
    : logic_(logic)
{
    CELER_EXPECT(!logic_.empty());
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
