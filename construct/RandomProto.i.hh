//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RandomProto.i.hh
 * \brief RandomProto inline method definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * True if fully defined
 */
RandomProto::Particle::operator bool() const
{
    return proto && volume_fraction > 0 && md;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
