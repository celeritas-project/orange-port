//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/ShapeTest.hh
 * \brief ShapeTest class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Nemesis/gtest/Test.hh"

#include <cmath>
#include <string>
#include <vector>
#include "orange/BoundingBox.hh"
#include "orange/Transform.hh"
#include "orange/construct/CSGCell.hh"
#include "orange/surfaces/SurfaceContainer.hh"

namespace celeritas
{
class Shape;

//---------------------------------------------------------------------------//
/*!
 * Unit test harness for ORANGE shapes.
 *
 * Its primary function is to create and execute a shape builder to do shape
 * comparisons. The typedefs and 'print_expected' are to aid testing.
 */
class ShapeTest : public ::Test
{
  public:
    //@{
    //! Public type aliases
    using Transform   = ::geometria::Transform;
    using Real3       = ::geometria::Real3;
    using SurfaceId   = ::celeritas::SurfaceId;
    using SurfaceType = ::celeritas::SurfaceType;
    using CSGCell     = ::celeritas::CSGCell;
    using VecDbl      = std::vector<real_type>;
    using VecStr      = std::vector<std::string>;
    using surfid_int  = SurfaceId::size_type;
    using logic_int   = CSGCell::logic_int;
    //@}

    enum LogicToken : CSGCell::logic_int
    {
        LOGIC_NOT = CSGCell::LOGIC_NOT,
        LOGIC_AND = CSGCell::LOGIC_AND,
        LOGIC_OR  = CSGCell::LOGIC_OR,
    };

    static constexpr real_type inf = HUGE_VAL;

  public:
    // Build a shape (with transformation)
    void build(const Shape& shape, Transform t);

    // Build a shape (no transformation)
    void build(const Shape& shape) { return this->build(shape, Transform{}); }

    // Print obtained values
    void print_expected() const;

    // Get the surface using an implicit value type, with error checking
    template<class S>
    S get_surface(surfid_int idx) const
    {
        EXPECT_LT(idx, this->surfaces.size());
        SurfaceId id(idx);
        EXPECT_EQ(this->surfaces.get_type(id), S::surface_type());
        return this->surfaces.get<S>(id);
    }

    // Get the surface data, checking for the correct surface type
    VecDbl get_surface_data(SurfaceType type, surfid_int idx) const;

    //// DATA ////

    SurfaceContainer surfaces; //!< Constructed surfaces
    CSGCell          cell;     //!< Region definition
    BoundingBox      bbox;     //!< Bounds in post-transform coordinate system
    VecStr surface_names;      //!< Surface names, in same order as in logic
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
