//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CSGCell.cc
 * \brief CSGCell class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CSGCell.hh"

#include <algorithm>
#include "base/FastHashSet.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Build a cell from a C string
 *
 * This is mostly used in unit tests.. A valid string satisfies the regex
 * "[0-9~!| ]+", but the result may not be a valid logic expression. The
 * constructor will ensure that the logic expression at least is consistent for
 * a CSG region definition (final 'stack depth' of 1)
 *
 * Example:
 * \code

     vols.insert(build_logic("4 ~ 5 & 6 &"), attr);

   \endcode
 */
CSGCell CSGCell::from_string(const char* c)
{
    VecLogic  logic;
    logic_int s = 0;
    while (char v = *c++)
    {
        if (v >= '0' && v <= '9')
        {
            // Parse a surface number. 'Push' this digit onto the surface ID by
            // multiplying the existing ID by 10.
            s = 10 * s + (v - '0');

            const char next = *c;
            if (next == ' ' || next == '\0')
            {
                // Next char is end of word or end of string
                logic.push_back(s);
                s = 0;
            }
        }
        else
        {
            // Parse a logic token
            switch (v)
            {
                // clang-format off
                case ' ': break;
                case '*': logic.push_back(LOGIC_TRUE); break;
                case '|': logic.push_back(LOGIC_OR);   break;
                case '&': logic.push_back(LOGIC_AND);  break;
                case '~': logic.push_back(LOGIC_NOT);  break;
                default:   CELER_ASSERT_UNREACHABLE();
                    // clang-format on
            }
        }
    }

    return CSGCell(std::move(logic));
}

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Construct from logic
 */
CSGCell::CSGCell(VecLogic logic) : logic_(std::move(logic))
{
    CELER_EXPECT(!logic_.empty());

    // Calculate max depth
    int max_depth = 0;
    int cur_depth = 0;

    for (auto id : logic_)
    {
        if (!is_operator_token(id) || id == LOGIC_TRUE)
        {
            ++cur_depth;
        }
        else if (id == LOGIC_AND || id == LOGIC_OR)
        {
            max_depth = std::max(cur_depth, max_depth);
            --cur_depth;
        }
    }
    logic_depth_ = std::max(cur_depth, max_depth);
    CELER_ASSERT(cur_depth == 1);

    CELER_ENSURE(logic_depth_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * Get sorted surface IDs used by this region definition
 */
auto CSGCell::get_local_surfids() const -> VecSurfaceId
{
    // Allocate set: number of operators has to be at least half of the logic
    // definition
    FastHashSet<SurfaceId> ids;
    ids.reserve((logic_.size() + 1) / 2);

    // Append items of the logic vector that aren't operator tokens
    for (logic_int surfid : logic_)
    {
        if (!is_operator_token(surfid))
            ids.insert(SurfaceId(surfid));
    }

    // Convert to unique sorted vector
    VecSurfaceId result(ids.begin(), ids.end());
    std::sort(result.begin(), result.end());
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Get the printable character corresponding to a logic token.
 */
char to_char(CSGCell::LogicToken t)
{
    CELER_EXPECT(t >= CSGCell::LOGIC_TRUE && t <= CSGCell::LOGIC_NOT);
    static const char symbols[] = {'*', '|', '&', '~'};
    return symbols[t - CSGCell::LOGIC_TRUE];
}

//---------------------------------------------------------------------------//
/*!
 * Write to stream
 */
std::ostream& operator<<(std::ostream& os, const CSGCell& cell)
{
    CELER_EXPECT(cell);

    for (CSGCell::logic_int t : cell.logic())
    {
        if (CSGCell::is_operator_token(t))
        {
            os << to_char(static_cast<CSGCell::LogicToken>(t));
        }
        else
        {
            os << t; // Surface ID integer
        }
        os << ' ';
    }
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
