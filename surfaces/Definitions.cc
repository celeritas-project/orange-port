//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/Definitions.cc
 * \brief SurfaceType definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Definitions.hh"
#include "base/Assert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Convert surface sense to equivalent C string
 */
const char* to_cstring(SignedSense type)
{
#define ORANGE_SIGNED_SENSE_CASE(VAL) \
    case SignedSense::VAL:            \
        return #VAL

    switch (type)
    {
        ORANGE_SIGNED_SENSE_CASE(inside);
        ORANGE_SIGNED_SENSE_CASE(on);
        ORANGE_SIGNED_SENSE_CASE(outside);
    }

#undef ORANGE_SIGNED_SENSE_CASE
    CELER_ASSERT_UNREACHABLE();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert surface type to equivalent C string
 */
const char* to_cstring(SurfaceType type)
{
    static const char* const names[] = {"px",
                                        "py",
                                        "pz",
                                        "so",
                                        "cxo",
                                        "cyo",
                                        "czo",
                                        "p",
                                        "s",
                                        "cx",
                                        "cy",
                                        "cz",
                                        "kx",
                                        "ky",
                                        "kz",
                                        "sq",
                                        "gq"};
    static_assert(
        sizeof(names) == sizeof(char*) * static_cast<int>(SurfaceType::size_),
        "SurfaceType string size mismatch");
    const char* name = names[static_cast<int>(type)];
    CELER_ENSURE(name);
    return name;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
