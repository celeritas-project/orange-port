//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CylCentered.pt.cc
 * \brief CylCentered explicit instantiations
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CylCentered.t.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//

template class CylCentered<Axis::x>;
template class CylCentered<Axis::y>;
template class CylCentered<Axis::z>;

template std::ostream& operator<<(std::ostream&, const CylCentered<Axis::x>&);
template std::ostream& operator<<(std::ostream&, const CylCentered<Axis::y>&);
template std::ostream& operator<<(std::ostream&, const CylCentered<Axis::z>&);

//---------------------------------------------------------------------------//
} // namespace celeritas
