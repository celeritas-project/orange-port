//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/SenseContainer.cc
 * \brief SenseContainer class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SenseContainer.hh"

#include <iostream>
#include "base/Join.hh"
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Print a sense container.
 */
std::ostream& operator<<(std::ostream& os, const SenseContainer& senses)
{
    os << '{' << join(senses.begin(), senses.end(), ' ', [](bool s) {
        return to_char(to_sense(s));
    }) << '}';
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
