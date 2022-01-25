//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/ObjectMetadata.hh
 * \brief ObjectMetadata class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <string>
#include "base/Provenance.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Simple metadata storage for an object.
 */
class ObjectMetadata
{
  public:
    //@{
    //! Public type aliases
    using SPConstProvenance = std::shared_ptr<const Provenance>;
    //@}

    struct Params
    {
        std::string       name;        //!< Name
        std::string       description; //!< Optional longer description
        SPConstProvenance provenance;  //!< Where the object was defined
    };

  public:
    // Construct default name/provenance from unknown object type
    static ObjectMetadata from_untitled_instance_of(const char* clsname);

    // Construct from parameters
    explicit ObjectMetadata(Params params);

    // Unassigned constructor
    ObjectMetadata() = default;

    //! Descriptive name (never empty)
    const std::string& name() const
    {
        CELER_EXPECT(*this);
        return data_.name;
    }

    //! Long description (maybe empty)
    const std::string& description() const { return data_.description; }

    // Provenance (originating file/line; always set)
    inline const SPConstProvenance& provenance() const;

    //! Whether metadata has been assigned
    explicit operator bool() const { return !data_.name.empty(); }

  private:
    //// DATA ////

    Params data_;
};

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
std::ostream& operator<<(std::ostream& os, const ObjectMetadata& md);

//---------------------------------------------------------------------------//
} // namespace celeritas

//! Construct object metadata for testing
#define ORANGE_MD_FROM_SOURCE(NAME) \
    (::celeritas::ObjectMetadata({NAME, {}, NEMESIS_PROVENANCE_FROM_SOURCE()}))

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "ObjectMetadata.i.hh"
//---------------------------------------------------------------------------//
