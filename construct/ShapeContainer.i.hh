//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ShapeContainer.i.hh
 * \brief ShapeContainer inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Assert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Add a shape, return reference to inserted shape
 */
auto ShapeContainer::insert(mapped_type shape) -> const mapped_type&
{
    CELER_EXPECT(shape);

    auto               insertion = shapes_.emplace(shape->name(), shape);
    const mapped_type& result    = insertion.first->second;
    Insist(insertion.second,
           "Failed to store duplicate shape '"
               << shape->name() << "' at " << *(shape->metadata().provenance())
               << ": previous shape at " << *(result->metadata().provenance()));

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Get the shape with the given name
 */
auto ShapeContainer::find(const key_type& name) const -> const_iterator
{
    return shapes_.find(name);
}

//---------------------------------------------------------------------------//
/*!
 * Get the shape with the given name
 */
auto ShapeContainer::at(const key_type& name) const -> const mapped_type&
{
    const_iterator found = this->find(name);
    Insist(found != shapes_.end(), "No shape named '" << name << "' is stored");
    return found->second;
}

//---------------------------------------------------------------------------//
/*!
 * Get the shape with the given name
 */
auto ShapeContainer::operator[](const key_type& name) const
    -> const mapped_type&
{
    return this->at(name);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
