//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RandomProto.hh
 * \brief RandomProto class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Proto.hh"

#include <string>
#include <vector>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Proto-universe for building universes packed with spheres.
 *
 * This simple proto-universe has only a single shape definition (external
 * boundary shape).
 */
class RandomProto : public Proto
{
  public:
    //@{
    //! Public type aliases
    using matid_type = geometria::matid_type;
    //@}

    struct Particle
    {
        SPConstProto   proto;
        real_type      volume_fraction{};
        ObjectMetadata md;

        // True if fully defined
        explicit inline operator bool() const;
    };

    struct Options
    {
        int seed                    = 0;    //!< Seed is hashed with proto name
        int failure_batch_size      = 1024; //!< Attempts for failure detection
        real_type failure_tolerance = 0.9;  //!< Fractional failures before
                                            //!< abort
    };

    struct Params
    {
        RegionVec             interior;
        matid_type            fill_matid{geometria::invalid_matid()};
        std::vector<Particle> particles;
        Options               options;
        ObjectMetadata        md;

        // True if fully defined
        explicit inline operator bool() const;
    };

  public:
    // Constructor
    explicit RandomProto(Params params);

    // Default destructor
    ~RandomProto();

    //! Options database
    const Options& options() const { return data_.options; }

    //! Get the proto metadata
    const ObjectMetadata& metadata() const final { return data_.md; }

    //! Get the boundary definition for defining a hole in a higher level
    const RegionVec& interior() const final { return data_.interior; }

    // Construct a universe from this proto
    BuildResult build(BuildArgs args) const final;

  private:
    //// DATA ////

    Params data_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

#include "RandomProto.i.hh"

//---------------------------------------------------------------------------//
