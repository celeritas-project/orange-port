//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file Fuzziness.hh
 * \brief Fuzziness class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Manage numerical soft equivalence used in geometry construction.
 *
 * There are tolerances that need to be set in construction for surface
 * deduplication or plane construction. We also have DBC checks that do soft
 * equivalence. This class is to make sure they're generally consistent with
 * each other.
 *
 * We generally have one or both of two kinds of "fuzziness". "Relative"
 * compares or expands compared to the relative size of an object or a particle
 * position; "absolute" compares to a fixed value.
 *
 * Currently all fuzziness parameters are set at the same time, using the
 * "length scale" (which should be the characteristic feature size of the
 * problem; e.g. ~1cm for a reactor or ~1e-6 cm if doing cellular irradiation
 * experiments). A global relative tolerance is given by the user and applied
 * to all relative tolerances, and an absolute tolerance is calculated by
 * multiplying the relative error into the length scale.
 *
 * \note Use the \c fuzziness() free function for ORANGE-specific units.
 * \todo: replace GG Fuzziness with an instance of this class.
 */
class Fuzziness
{
  public:
    // Construct with tolerance and length scale
    Fuzziness(real_type tolerance, real_type length_scale);

    //! Construct with default length scale of 1 cm
    explicit Fuzziness(real_type eps) : Fuzziness(eps, 1.0) {}

    //! Construct with default tolerance and length scale of 1 cm
    Fuzziness() : Fuzziness(1e-8, 1.0) {}

    //// GLOBAL ASSIGNMENT ////

    //! Get the characteristic problem length scale
    real_type length_scale() const { return length_scale_; }

    //// SHAPES/SURFACE CONSTRUCTION ////

    //! Relative error for surface elision
    real_type surface_elision_rel() const { return surface_elision_rel_; }
    //! Absolute error for surface elision
    real_type surface_elision_abs() const { return surface_elision_abs_; }
    //! Absolute error for dropping higher-order quadric coefficients
    real_type surface_simplification_abs() const
    {
        return surface_elision_abs_;
    }

    //! Relative error when testing shape overlaps
    real_type shape_enclosure_rel() const { return shape_enclosure_rel_; }

    //! Bounding box tolerance (expand by this fraction)
    real_type bbox_rel() const { return bbox_rel_; }

    //// TRANSPORT ////

    //! Cell initialization bump distance
    real_type bump_rel() const { return bump_rel_; }
    //! Absolute bump distance, done in cell initialization
    real_type bump_abs() const { return bump_abs_; }

    //! Minimum (slightly negative) distance allowed during transport
    real_type min_dist() const { return min_dist_; }

    //// SURFACE CALCULATIONS ////

    //! Tolerance for directions being along a surface
    static constexpr real_type quadratic_abs() { return 1.e-10; }

  private:
    // Characteristic length scale of the problem in cm
    real_type length_scale_;

    // Relative and absolute values for surface elision
    real_type surface_elision_abs_;
    real_type surface_elision_rel_;

    // Relative enclosure test for shape intersection
    real_type shape_enclosure_rel_;

    // Bounding boxes are expanded by this much during BVH construction
    real_type bbox_rel_;

    // Relative and absolute initialization 'bump' distances
    real_type bump_rel_;
    real_type bump_abs_;

    // Minimum (slight negative) distance allowed during transport
    real_type min_dist_;
};

//---------------------------------------------------------------------------//
/*!
 * Instance of fuzziness used by the ORANGE package.
 *
 * Ideally this is be owned by a geometry instance and shared by classes
 * that need it: tracking errors might result from constructing with one
 * tolerance and transporting on another.
 */
Fuzziness& fuzziness();

//---------------------------------------------------------------------------//
} // namespace celeritas
