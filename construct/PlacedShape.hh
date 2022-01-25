//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/PlacedShape.hh
 * \brief PlacedShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <memory>
#include <string>
#include "base/Provenance.hh"
#include "orange/Transform.hh"
#include "../Definitions.hh"
#include "orange/query/ObjectMetadata.hh"
#include "Shape.hh"

namespace celeritas
{
class UnitBuilder;

//---------------------------------------------------------------------------//
/*!
 * Shape with a transformation and metadata.
 *
 * Shapes are simple classes used to create surfaces. The "placed" shape allows
 * these surfaces to be defined from a transformed coordinate system, and it
 * allows the shape to be tagged with other metadata like a name and
 * provenance.
 */
class PlacedShape
{
  public:
    //@{
    //! Public type aliases
    using SPConstShape      = std::shared_ptr<const Shape>;
    using SPConstProvenance = std::shared_ptr<const Provenance>;
    //@}

    //! User-facing properties of this shape
    struct Params
    {
        SPConstShape   shape;     //!< Core shape definition
        Transform      transform; //!< Applied transformation
        ObjectMetadata md;        //!< Name and provenance
    };

  public:
    // Construct with fixed parameters
    explicit PlacedShape(Params params);

    //// ACCESSORS ////

    //! Shape class (always valid)
    const SPConstShape& shape() const { return data_.shape; }

    //! Transformation
    const Transform& transform() const { return data_.transform; }

    //! Descriptive name (never empty)
    const std::string& name() const { return data_.md.name(); }

    //! All metadata
    const ObjectMetadata& metadata() const { return data_.md; }

    //// CONSTRUCTION ////

    // Return a transformed version of this shape
    PlacedShape transformed(const Transform& t) const;

  private:
    Params data_;
};

//---------------------------------------------------------------------------//
// Print a shape
std::ostream& operator<<(std::ostream&, const PlacedShape&);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
