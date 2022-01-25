//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/Shape.cc
 * \brief Shape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Shape.hh"

#include <iostream>

namespace celeritas
{
//---------------------------------------------------------------------------//
//! Default virtual destructor
Shape::~Shape() = default;

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Print a shape.
 *
 * \todo Each shape has additional information
 */
std::ostream& operator<<(std::ostream& os, const Shape& shape)
{
    os << '{' << shape.type() << '}';
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
