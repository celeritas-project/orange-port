//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/UnitProto.i.hh
 * \brief UnitProto inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * True if fully defined
 */
UnitProto::Media::operator bool() const
{
    return !interior.empty() && matid != geometria::invalid_matid() && md;
}

//---------------------------------------------------------------------------//
/*!
 * True if fully defined
 */
UnitProto::Array::operator bool() const
{
    return proto && !interior.empty() && md;
}

//---------------------------------------------------------------------------//
/*!
 * True if fully defined
 */
UnitProto::Hole::operator bool() const
{
    return proto && md;
}

//---------------------------------------------------------------------------//
/*!
 * \brief True if fully defined
 */
UnitProto::Boundary::operator bool() const
{
    return !interior.empty() && md;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
