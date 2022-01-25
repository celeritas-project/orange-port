//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/FaceNameCalculator.hh
 * \brief FaceNameCalculator class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <vector>

namespace celeritas
{
class Shape;
//---------------------------------------------------------------------------//
/*!
 * Get the (ordered) names of the surfaces of a shape.
 */
class FaceNameCalculator
{
  public:
    //@{
    //! Public type aliases
    using result_type = std::vector<std::string>;
    //@}

  public:
    // Determine face names
    result_type operator()(const Shape& shape);
};

//---------------------------------------------------------------------------//
} // namespace celeritas
