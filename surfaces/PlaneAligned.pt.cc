//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/PlaneAligned.pt.cc
 * \brief PlaneAligned explicit instantiations
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "PlaneAligned.t.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//

template class PlaneAligned<Axis::x>;
template class PlaneAligned<Axis::y>;
template class PlaneAligned<Axis::z>;

template std::ostream& operator<<(std::ostream&, const PlaneAligned<Axis::x>&);
template std::ostream& operator<<(std::ostream&, const PlaneAligned<Axis::y>&);
template std::ostream& operator<<(std::ostream&, const PlaneAligned<Axis::z>&);

//---------------------------------------------------------------------------//
} // namespace celeritas
