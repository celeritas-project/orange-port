//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ShapeContainer.cc
 * \brief ShapeContainer class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ShapeContainer.hh"

#include <iostream>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Print to a stream for debugging.
 */
std::ostream& operator<<(std::ostream& os, const ShapeContainer& shapes)
{
    os << '{';
    for (auto iter = shapes.cbegin(); iter != shapes.cend(); ++iter)
    {
        os << *iter->second << ",\n";
    }
    os << '}';
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
