//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file Fuzziness.cc
 * \brief Fuzziness class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fuzziness.hh"

#include "base/Assert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with tolerance and length scale.
 *
 * Sets all relative tolerances to the given eps, absolute tolerances to
 * eps length_scale, and the minimum bump distance to -eps.
 *
 * The length scale of the problem is its characteristic length [cm].
 */
Fuzziness::Fuzziness(real_type eps, real_type length_scale = 1.0)
{
    CELER_VALIDATE(eps > 0 && eps < 1,
                   << "Tolerance must be between 0 and 1; given " << eps);
    CELER_VALIDATE(length_scale > 0,
                   << "Length scale must be positive: given " << length_scale);

    length_scale_ = length_scale;

    surface_elision_abs_ = eps * length_scale_;
    surface_elision_rel_ = eps;
    shape_enclosure_rel_ = eps;
    bbox_rel_            = eps * 10.0;
    bump_rel_            = eps;
    bump_abs_            = eps * length_scale_;

    // Minimum allowed distance-to-surface gets a slightly negative value
    min_dist_ = -eps * length_scale_;

    CELER_ENSURE(min_dist_ < 0);
}

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Instance of fuzziness used by the ORANGE package.
 *
 * The fuzziness can be changed by assigning a new one:
 * \code
 * celeritas::fuzziness() = celeritas::Fuzziness{1e-5};
 * \endcode
 */
Fuzziness& fuzziness()
{
    static Fuzziness fuzz{};
    return fuzz;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
