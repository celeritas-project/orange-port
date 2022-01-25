//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/SenseContainer.i.hh
 * \brief SenseContainer inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Macros.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
SenseContainer::SenseContainer(std::initializer_list<value_type> values)
    : storage_(values)
{
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::operator[](size_type index) -> reference
{
    return storage_[index];
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::operator[](size_type index) const -> const_reference
{
    return storage_[index];
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::begin() -> iterator
{
    return storage_.begin();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::end() -> iterator
{
    return storage_.end();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::begin() const -> const_iterator
{
    return storage_.begin();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::end() const -> const_iterator
{
    return storage_.end();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::cbegin() const -> const_iterator
{
    return storage_.cbegin();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::cend() const -> const_iterator
{
    return storage_.cend();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
bool SenseContainer::empty() const
{
    return storage_.empty();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
auto SenseContainer::size() const -> size_type
{
    return storage_.size();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
void SenseContainer::clear()
{
    return storage_.clear();
}

//---------------------------------------------------------------------------//
CELER_FORCEINLINE_FUNCTION
void SenseContainer::resize(size_type size)
{
    return storage_.resize(size);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
