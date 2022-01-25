//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CSGCell.hh
 * \brief CSGCell class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <type_traits>
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Constructive solid geometry definition of a region by surfaces.
 *
 * Currently, CSG cells are represented as a flattened Reverse Polish
 * Notation (RPN) expression where each input value is a surface ID or a
 * special "logic token".
 *
 * This replaces the \c gg::Volume class.
 */
class CSGCell
{
  public:
    //@{
    //! Public type aliases
    using logic_int    = std::make_signed<SurfaceId::size_type>::type;
    using VecLogic     = std::vector<logic_int>;
    using VecSurfaceId = std::vector<SurfaceId>;
    //@}

    //! Logic relating surface definitions.
    // Negative integers are reserved for special tokens.
    enum LogicToken : logic_int
    {
        BEGIN_LOGIC_TOKEN = logic_int(-4),
        LOGIC_TRUE        = BEGIN_LOGIC_TOKEN, //!< Invariant 'true'
        LOGIC_OR,                              //!< Binary logical OR
        LOGIC_AND,                             //!< Binary logical AND
        LOGIC_NOT,                             //!< Unary negation
    };

    //! Whether an integer is a special logic token
    static bool is_operator_token(logic_int lv) { return (lv < 0); }

  public:
    // Construct from RPN string, used by unit tests
    static CSGCell from_string(const char* logic);

    // Default empty constructor
    CSGCell() = default;

    // Construct from logic
    CSGCell(VecLogic logic);

    // Whether this cell has been constructed
    explicit operator bool() const { return !logic_.empty(); }

    //// ACCESSORS ////

    //! Get the surface logic representation
    const VecLogic& logic() const
    {
        CELER_EXPECT(*this);
        return logic_;
    }

    // Get sorted surface IDs used by this CSG voume
    VecSurfaceId get_local_surfids() const;

    //! Levels in the CSG tree / max logic stack depth
    int logic_depth() const { return logic_depth_; }

  private:
    VecLogic logic_;
    int      logic_depth_ = 0;
};

//---------------------------------------------------------------------------//
// Get the printable character corresponding to a logic token
char to_char(CSGCell::LogicToken);

//---------------------------------------------------------------------------//
// Write to stream
std::ostream& operator<<(std::ostream& os, const CSGCell& cell);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
