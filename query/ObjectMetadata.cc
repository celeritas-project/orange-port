//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/ObjectMetadata.cc
 * \brief ObjectMetadata class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ObjectMetadata.hh"

#include <iostream>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct default name/provenance from unknown object type
 */
ObjectMetadata ObjectMetadata::from_untitled_instance_of(const char* clsname)
{
    Params params;
    params.name       = std::string("Untitled ") + clsname;
    params.provenance = NEMESIS_PROVENANCE_FROM_SOURCE_DEBUG();

    return ObjectMetadata(std::move(params));
}

//---------------------------------------------------------------------------//
/*!
 * Construct from parameters
 */
ObjectMetadata::ObjectMetadata(Params params) : data_(std::move(params))
{
    CELER_EXPECT(!data_.name.empty());
    CELER_EXPECT(data_.provenance);
    CELER_ENSURE(*this);
}

//---------------------------------------------------------------------------//
/*!
 * Print metadata
 */
std::ostream& operator<<(std::ostream& os, const ObjectMetadata& md)
{
    CELER_EXPECT(md);

    os << '"' << md.name() << '"';
    if (!md.description().empty())
    {
        os << ": " << md.description();
    }
    os << " (" << *md.provenance() << ")";
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
