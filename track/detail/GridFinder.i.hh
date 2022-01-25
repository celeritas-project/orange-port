//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/GridFinder.i.hh
 * \brief GridFinder inline method definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <algorithm>
#include "base/Assert.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Construct with a sorted grid.
 */
template<class T>
GridFinder<T>::GridFinder(span<const T> grid) : grid_(grid)
{
    CELER_EXPECT(grid_.size() >= 2);
    CELER_EXPECT(grid_.front() <= grid_.back()); // "sorted"
}

//---------------------------------------------------------------------------//
/*!
 * Find on or around the grid.
 */
template<class T>
auto GridFinder<T>::operator()(argument_type val) const -> result_type
{
    result_type result;

    auto iter = std::lower_bound(grid_.begin(), grid_.end(), val);
    if (iter == grid_.begin())
    {
        // Below first point
        result.edge = -1;
    }
    else if (iter == grid_.end() || (iter == grid_.end() - 1 && val == *iter))
    {
        // Above or on final grid point: move to high end of last bin
        iter        = grid_.end() - 2;
        result.edge = 1;
    }
    else if (val == *iter)
    {
        // On a grid point
        result.edge = -1;
    }
    else
    {
        // Between two grid points
        --iter;
        result.edge = 0;
    }
    result.cell = iter - grid_.begin();

    CELER_ENSURE(result.edge >= -1 && result.edge <= 1);
    CELER_ENSURE(result.cell + 1 < grid_.size());
    return result;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
