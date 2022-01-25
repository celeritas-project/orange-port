//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/SurfaceInserter.i.hh
 * \brief SurfaceInserter inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Add a surface (type-deleted) with the given coefficients
 */
template<class S>
SurfaceId SurfaceInserter::operator()(const S& surface)
{
    return (*this)(surface.view());
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
