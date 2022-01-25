//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CylAligned.pt.cc
 * \brief CylAligned explicit instantiations
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CylAligned.t.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//

template class CylAligned<Axis::x>;
template class CylAligned<Axis::y>;
template class CylAligned<Axis::z>;

template std::ostream& operator<<(std::ostream&, const CylAligned<Axis::x>&);
template std::ostream& operator<<(std::ostream&, const CylAligned<Axis::y>&);
template std::ostream& operator<<(std::ostream&, const CylAligned<Axis::z>&);

//---------------------------------------------------------------------------//
} // namespace celeritas
