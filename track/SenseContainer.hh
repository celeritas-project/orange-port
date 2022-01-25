//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/SenseContainer.hh
 * \brief SenseContainer class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <vector>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Store surface sense state inside a cell.
 *
 * This is currently a subset of std::vector functionality.
 *
 * \todo Change value_type to Sense
 */
class SenseContainer
{
  public:
    //@{
    //! Public type aliases
    using Storage_t       = std::vector<bool>;
    using value_type      = bool;
    using size_type       = Storage_t::size_type;
    using iterator        = Storage_t::iterator;
    using const_iterator  = Storage_t::const_iterator;
    using reference       = Storage_t::reference;
    using const_reference = Storage_t::const_reference;
    //@}

  public:
    //@{
    //! Constructors
    SenseContainer() = default;
    inline SenseContainer(std::initializer_list<value_type>);
    //@}

    //@{
    //! Element access
    inline reference       operator[](size_type);
    inline const_reference operator[](size_type) const;
    //@}

    //@{
    //! Iterators
    inline iterator       begin();
    inline iterator       end();
    inline const_iterator begin() const;
    inline const_iterator end() const;
    inline const_iterator cbegin() const;
    inline const_iterator cend() const;
    //@}

    //@{
    //! Capacity
    inline bool      empty() const;
    inline size_type size() const;
    //@}

    //@{
    //! Modifiers
    inline void clear();
    inline void resize(size_type);
    //@}

  private:
    //// DATA ////
    std::vector<bool> storage_;
};

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//

std::ostream& operator<<(std::ostream& os, const SenseContainer& senses);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "SenseContainer.i.hh"
//---------------------------------------------------------------------------//
