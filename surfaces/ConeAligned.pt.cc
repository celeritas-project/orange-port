//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/ConeAligned.pt.cc
 * \brief ConeAligned explicit instantiations
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ConeAligned.t.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//

template class ConeAligned<Axis::x>;
template class ConeAligned<Axis::y>;
template class ConeAligned<Axis::z>;

template std::ostream& operator<<(std::ostream&, const ConeAligned<Axis::x>&);
template std::ostream& operator<<(std::ostream&, const ConeAligned<Axis::y>&);
template std::ostream& operator<<(std::ostream&, const ConeAligned<Axis::z>&);

//---------------------------------------------------------------------------//
} // namespace celeritas
