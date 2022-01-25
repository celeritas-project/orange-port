//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/LogicStack.cc
 * \brief LogicStack class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "LogicStack.hh"

#include <iostream>

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Print logic stack
 */
std::ostream& operator<<(std::ostream& os, const LogicStack& stack)
{
    os << '{';
    for (int i = 0; i != stack.size(); ++i)
    {
        os << (stack[i] ? 'T' : 'F') << ' ';
    }
    os << '}';
    return os;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
