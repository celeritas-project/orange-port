//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/LogicStack.hh
 * \brief LogicStack class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <cstdint>
#include <iosfwd>

#include "base/Macros.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Simple, fixed-max-size (no heap allocation) stack.
 *
 * This uses a bit field where the "top" of the stack is the least significant
 * bit.
 *
 * For a KENO-type region definition, the max stack depth is 2. This logic
 * stack can hold up to 32 entries (minimum size of storage_type).
 *
 * The underlying code is highly optimized. For example, calling 'apply_and'
 * inside the 'evaluate' function (on GCC 5 with -O2) results in: \verbatim
        movl    %eax, %ecx  # stack, temp
        shrl    %eax        # D.44391
        andl    $1,   %ecx  #, temp
        andl    %ecx, %eax  # temp, stack
 * \endverbatim
 *
 * Furthermore, and delightfully, if LogicStack is local to a function, and
 * DBC is off, all operations on size_ are optimized out completely! So
 * there's no penalty to adding that extra safety check.
 */
class LogicStack
{
  public:
    //@{
    //! Typedefs
    using value_type   = bool;
    using storage_type = std::uint_fast32_t;
    using size_type    = storage_type;
    //@}

  public:
    CELER_FORCEINLINE_FUNCTION LogicStack() = default;

    //// ACCESSORS ////

    //! Number of elements on the stack
    CELER_FORCEINLINE_FUNCTION size_type size() const { return size_; }

    // Whether any elements exist
    inline bool empty() const;

    // Access the top value of the stack
    inline value_type top() const;

    // Access a single bit (zero is deepest level of stack), used by ostream
    inline value_type operator[](size_type index) const;

    //// MUTATORS ////

    // Push a boolean onto the stack
    inline void push(value_type v);

    // Pop a value off the stack
    inline value_type pop();

    // Negate the value on the top of the stack
    inline void apply_not();

    // Apply boolean 'and' to the top of the stack
    inline void apply_and();

    // Apply boolean 'or' to the top of the stack
    inline void apply_or();

  private:
    //! Get the least significant bit
    CELER_FORCEINLINE_FUNCTION static constexpr storage_type
    lsb(storage_type val)
    {
        return val & storage_type(1);
    }

    //! Shift right by one
    CELER_FORCEINLINE_FUNCTION static constexpr storage_type
    shr(storage_type val)
    {
        return val >> storage_type(1);
    }

    //! Shift left by one
    CELER_FORCEINLINE_FUNCTION static constexpr storage_type
    shl(storage_type val)
    {
        return val << storage_type(1);
    }

    //! Greatest number of boolean values allowed on the stack
    CELER_FORCEINLINE_FUNCTION static constexpr storage_type max_stack_depth()
    {
        return sizeof(storage_type) * 8;
    };

  private:
    //// DATA ////

    // Stack data
    storage_type data_ = 0;

    // Depth
    storage_type size_ = 0;
};

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//

// Print logic stack
std::ostream& operator<<(std::ostream& os, const LogicStack& stack);

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "LogicStack.i.hh"
//---------------------------------------------------------------------------//
