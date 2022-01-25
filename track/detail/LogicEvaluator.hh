//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/LogicEvaluator.hh
 * \brief LogicEvaluator class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Nemesis/containers/Span.hh"
#include "orange/Definitions.hh"
#include "../SenseContainer.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Evaluate a logical expression applied to a vector of bools.
 */
class LogicEvaluator
{
  public:
    //@{
    //! Public type aliases
    using SpanConstLogic = span<const logic_int>;
    //@}

    //! Special logical Evaluator tokens.
    // The enum values are set to the highest 4 values of logic_int.
    enum LogicToken : logic_int
    {
        logic_begin = logic_int(~logic_int(4)),
        logic_true  = logic_begin, //!< Push 'true'
        logic_or,                  //!< Binary logical OR
        logic_and,                 //!< Binary logical AND
        logic_not,                 //!< Unary negation
        logic_end
    };

    //! Whether an integer is a special logic token
    static bool is_operator_token(logic_int lv)
    {
        return (lv >= logic_begin);
    }

    //! Largest number of surfaces used in a single logic def
    constexpr static size_type max_num_surfaces()
    {
        return static_cast<size_type>(logic_begin);
    }

  public:
    // Construct with a view to some logic definition
    explicit inline LogicEvaluator(SpanConstLogic logic);

    // Evaluate a logical expression, substituting bools from the vector
    bool operator()(const SenseContainer& values) const;

  private:
    //// DATA ////

    SpanConstLogic logic_;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "LogicEvaluator.i.hh"
//---------------------------------------------------------------------------//
