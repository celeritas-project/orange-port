//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/BBoxCalculator.hh
 * \brief BBoxCalculator class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "orange/BoundingBox.hh"
#include "orange/Transform.hh"

namespace celeritas
{
class PlacedShape;
class Shape;

//---------------------------------------------------------------------------//
/*!
 * Calculate bounding boxes for shapes.
 *
 * This is a simple helper class that can be used in downstream applications.
 */
class BBoxCalculator
{
  public:
    //@{
    //! Public type aliases
    using result_type = geometria::BoundingBox;
    using Transform   = geometria::Transform;
    //@}

  public:
    // Constructor
    BBoxCalculator() = default;

    // Calculate the bounding box for a transformed shape
    result_type operator()(const Shape&, const Transform& t);

    // Calculate the bounding box for a shape
    result_type operator()(const Shape&);

    // Calculate the bounding box for a placed shape
    result_type operator()(const PlacedShape&);
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
