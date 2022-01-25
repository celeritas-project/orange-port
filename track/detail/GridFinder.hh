//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/GridFinder.hh
 * \brief GridFinder class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Nemesis/containers/Span.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Find the cell and edge inside a grid.
 *
 * Examples:
 * \verbatim
     x----|----|   ==> {0,-1}
     |-x--|----|   ==> {0, 0}
     |----x----|   ==> {1,-1}
     |----|--x-|   ==> {1, 0}
     |----|----x   ==> {1, 1}
   x |----|----|   ==> {0,-1}
     |----|----|x  ==> {1, 1}
 * \endverbatim
 */
//---------------------------------------------------------------------------//

template<class T>
class GridFinder
{
  public:
    //!@{
    //! Public type aliases
    using argument_type = T;
    //!@}

    struct result_type
    {
        unsigned int cell; //!< [0, span.size())
        int          edge; //!< {-1, 0, or 1}

        bool operator==(const result_type& other) const
        {
            return cell == other.cell && edge == other.edge;
        }
    };

  public:
    // Construct with a sorted grid
    inline GridFinder(span<const T> grid);

    // Find on the grid
    inline result_type operator()(argument_type value) const;

  private:
    //// DATA ////
    span<const T> grid_;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "GridFinder.i.hh"
//---------------------------------------------------------------------------//
