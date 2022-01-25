//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/LogicEvaluator.cc
 * \brief LogicEvaluator class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "LogicEvaluator.hh"

#include "base/Assert.hh"
#include "LogicStack.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Evaluate a logical expression, substituting bools from the vector
 */
bool LogicEvaluator::operator()(const SenseContainer& values) const
{
    LogicStack stack;

    for (logic_int lgc : logic_)
    {
        if (lgc < logic_begin)
        {
            // Boolean from the value list; push onto the stack
            CELER_EXPECT(lgc < values.size());
            stack.push(values[lgc]);
            continue;
        }

        // Apply logic operator
        switch (lgc)
        {
            // clang-format off
            case logic_true: stack.push(true);  break;
            case logic_or:   stack.apply_or();  break;
            case logic_and:  stack.apply_and(); break;
            case logic_not:  stack.apply_not(); break;
            default:         CELER_ASSERT_UNREACHABLE();
                // clang-format on
        }
    }
    CELER_ENSURE(stack.size() == 1);
    return stack.top();
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
