//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/LogicStack.i.hh
 * \brief LogicStack inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Assert.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Whether any elements exist
 */
CELER_FORCEINLINE_FUNCTION bool LogicStack::empty() const
{
    return size_ == storage_type(0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access the top value of the stack
 */
CELER_FORCEINLINE_FUNCTION auto LogicStack::top() const -> value_type
{
    CELER_EXPECT(!empty());
    return LogicStack::lsb(data_);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a single bit (zero is deepest level of stack)
 */
auto LogicStack::operator[](size_type index) const -> value_type
{
    CELER_EXPECT(index < size());
    storage_type shift = size() - index - storage_type(1);
    return LogicStack::lsb(data_ >> shift);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Push a boolean onto the stack
 */
CELER_FORCEINLINE_FUNCTION void LogicStack::push(value_type v)
{
    CELER_EXPECT(size() != max_stack_depth());
    // Shift stack left and add least significant bit
    data_ = LogicStack::shl(data_) | LogicStack::lsb(v);
    // Size for DBC
    ++size_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pop a value off the stack
 */
CELER_FORCEINLINE_FUNCTION auto LogicStack::pop() -> value_type
{
    CELER_EXPECT(!empty());
    // Extract least significant bit
    value_type result = LogicStack::lsb(data_);
    // Shift right
    data_ = LogicStack::shr(data_);
    // Update size
    --size_;
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Negate the value on the top of the stack
 */
CELER_FORCEINLINE_FUNCTION void LogicStack::apply_not()
{
    CELER_EXPECT(!empty());
    data_ ^= storage_type(1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply boolean 'and' to the top of the stack
 */
CELER_FORCEINLINE_FUNCTION void LogicStack::apply_and()
{
    CELER_EXPECT(size() >= storage_type(2));
    storage_type temp = LogicStack::lsb(data_);
    data_             = LogicStack::shr(data_) & (temp | ~storage_type(1));
    --size_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply boolean 'or' to the top of the stack
 */
CELER_FORCEINLINE_FUNCTION void LogicStack::apply_or()
{
    CELER_EXPECT(size() >= storage_type(2));
    data_ = LogicStack::shr(data_) | LogicStack::lsb(data_);
    --size_;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
