//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ArrayProto.hh
 * \brief ArrayProto class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Proto.hh"

#include <cstddef>
#include "base/HyperslabVector.hh"
#include "base/Array.hh"
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Interface class for constructing a KENO array.
 */
class ArrayProto : public Proto
{
  public:
    //@{
    //! Public type aliases
    using DimVector = Array<size_type, 3>;
    using ProtoVec3 = HyperslabVector<SPConstProto, 3>;
    //@}

  public:
    // Transformation to move the 'low' corner to the given point
    virtual Transform calc_placement(SpanConstReal3 pos) const = 0;

    // Transformation to move the origin of the specified unit to the parent
    virtual Transform
    calc_placement(DimVector index, SpanConstReal3 pos) const = 0;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
