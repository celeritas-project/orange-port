//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/HashPair.hh
 * \brief HashPair class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <functional>
#include "base/HashFunctions.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * \struct HashPair
 * Create a well-distributed hash from a std::pair.
 */
template<class Pair_T>
struct HashPair
{
    using argument_type = const Pair_T&;
    using first_type    = typename Pair_T::first_type;
    using second_type   = typename Pair_T::second_type;
    using result_type   = size_type;

    result_type operator()(argument_type p) const
    {
        return hash_combine(hash_first(p.first), hash_second(p.second));
    }

    std::hash<first_type>  hash_first;
    std::hash<second_type> hash_second;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
