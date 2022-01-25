//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/TrackerTest.cc
 * \brief TrackerTest class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "TrackerTest.hh"

#include <algorithm>
#include <sstream>
#include "base/Range.hh"
#include "base/StringFunctions.hh"
#include "base/VectorFunctions.hh"
#include "orange/track/Tracker.hh"
#include "orange/query/UniverseMetadata.hh"
#include "../Tracker.hh"
#include "../TrackingError.hh"

using namespace celeritas;
using std::cout;
using std::endl;

namespace orange_test
{
#ifdef ORANGE_INTEL_CONSTEXPR_BUG
const real_type TrackerTest::inf = celeritas::no_intersection();
#endif

//---------------------------------------------------------------------------//
TrackerTest::TrackerTest()
{
    state_ref_.pos            = make_fixed_view(pos_);
    state_ref_.dir            = make_fixed_view(dir_);
    state_ref_.cell           = {};
    state_ref_.surface        = {};
    state_ref_.sense          = {};
    state_ref_.temp_senses    = &temp_senses_;
    state_ref_.temp_face_dist = &temp_face_dist_;
}

//---------------------------------------------------------------------------//
/*!
 * Initialize with a given pos/dir
 */
void TrackerTest::set_state(const Real3& pos,
                            const Real3& dir,
                            VolumeId     cell,
                            SurfaceId    surface,
                            Sense        sense)
{
    CELER_EXPECT(tracker);
    CELER_EXPECT(!cell || cell.get() < tracker->num_volumes());
    CELER_EXPECT(!surface || surface.get() < tracker->num_surfaces());

    pos_ = pos;
    dir_ = dir;
    normalize_direction(&dir_);
    state_ref_.cell    = cell;
    state_ref_.surface = surface;
    state_ref_.sense   = sense;

    CELER_ENSURE(state_ref_.temp_senses && state_ref_.temp_face_dist);
}

//---------------------------------------------------------------------------//
/*!
 * Find the cell from its label (nullptr allowed)
 */
VolumeId TrackerTest::find_cell(const char* label) const
{
    VolumeId volume_id;
    if (label)
    {
        auto iter = vol_ids_.find(label);
        CELER_VALIDATE(iter != vol_ids_.end(),
               << "nonexistent volume label '" << label << '\'');
        volume_id = iter->second;
    }
    return volume_id;
}

//---------------------------------------------------------------------------//
/*!
 * Find the surface from its label (nullptr allowed)
 */
SurfaceId TrackerTest::find_surface(const char* label) const
{
    SurfaceId surface_id;
    if (label)
    {
        auto iter = surf_ids_.find(label);
        CELER_VALIDATE(iter != surf_ids_.end(),
               << "nonexistent surface label '" << label << '\'');
        surface_id = iter->second;
    }
    return surface_id;
}

//---------------------------------------------------------------------------//
/*!
 * Cell name (or sentinel if no surface).
 */
std::string TrackerTest::id_to_label(VolumeId cell) const
{
    CELER_EXPECT(this->md);
    if (!cell)
        return "[none]";

    return this->md->id_to_label(cell);
}

//---------------------------------------------------------------------------//
/*!
 * Surface name (or sentinel if no surface).
 */
std::string TrackerTest::id_to_label(SurfaceId surf) const
{
    CELER_EXPECT(this->md);
    if (!surf)
        return "[none]";

    return this->md->id_to_label(surf);
}

//---------------------------------------------------------------------------//
/*!
 * Get the string output from the unit metadata.
 *
 * TODO when all our compilers support the regular expressions library, use
 * that instead.
 */
std::string TrackerTest::describe_md() const
{
    CELER_EXPECT(md && tracker);

    std::ostringstream os;
    os << '\n' << celeritas::to_stream(*md, *tracker);
    std::string result = os.str();

    // Basically equivalent to re.sub(r'``(.*?):(\d+)``', r'``\1:NNN``',
    // result)
    enum State
    {
        NORMAL,
        FIRST_BACKTICK,
        CODE,
        NUMBER,
        CODE_BACKTICK,
    } state;
    state = NORMAL;
    for (char& c : result)
    {
        switch (state)
        {
            case NORMAL:
                if (c == '`')
                    state = FIRST_BACKTICK;
                break;
            case FIRST_BACKTICK:
                if (c == '`')
                    state = CODE;
                else
                    state = NORMAL;
                break;
            case CODE:
                if (c == ':')
                    state = NUMBER;
                break;
            case NUMBER:
                if ('0' <= c && c <= '9')
                    c = 'N';
                else if (c == '`')
                    state = CODE_BACKTICK;
                break;
            case CODE_BACKTICK:
                if (c == '`')
                    state = NORMAL;
                break;
        };
    }

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Construct from tracker/MD
 */
void TrackerTest::set_tracker(UPConstTracker      tracker,
                              SPConstUnivMetadata univ_md)
{
    CELER_EXPECT(tracker);
    CELER_EXPECT(univ_md);
    this->tracker = std::move(tracker);
    this->md      = std::move(univ_md);

    const UniverseMetadata& md = *this->md;

    for (auto c : range(md.num_volumes()))
    {
        auto insertion = this->volume_ids_.insert(
            {md.id_to_label(VolumeId{c}), VolumeId{c}});
        CELER_ASSERT(insertion.second);
    }

    for (auto s : range(md.num_surfaces()))
    {
        auto insertion = surface_ids_.insert(
            {md.id_to_label(SurfaceId{s}), SurfaceId{s}});
        CELER_ASSERT(insertion.second);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Track until an infinite distance or exception is reached
 */
auto TrackerTest::track(const Real3& pos, const Real3& dir) -> TrackResult
{
    CELER_EXPECT(this->tracker && this->md);
    TrackResult result;

    // Initialize internal 'state'
    this->set_state(pos, dir);

    constexpr int max_zero_movements = 2;
    int           zero_movement      = 0;

    while (true)
    {
        Initialization init;
        try
        {
            init = tracker->initialize(this->state_ref());
        }
        catch (const TrackingError& e)
        {
            result.cells.push_back(e.what());
            if (auto* oce = dynamic_cast<const OverlappingCellError*>(&e))
            {
                for (auto cid : oce->local_vols())
                {
                    result.cells.push_back(this->md->id_to_label(cid));
                }
            }
            break;
        }
        catch (const std::exception& e)
        {
            ADD_FAILURE() << "Failed to initialize at " << pos_ << " along "
                          << dir_ << ": " << e.what() << endl;
            result.cells.push_back("[error]");
            break;
        }

        if (init.cell)
        {
            result.cells.push_back(this->md->id_to_label(init.cell));
        }
        else
        {
            // Failed to find cell; leave it to higher-level code to bump and
            // retry
            result.cells.push_back("[fail]");
            break;
        }

        if (init.surface)
        {
            result.surfaces.push_back(this->md->id_to_label(init.surface));
            result.senses.push_back(to_char(init.sense));
        }
        else
        {
            result.surfaces.push_back("");
            result.senses.push_back(' ');
        }

        if (state_ref_.cell == init.cell)
        {
            // Cell didn't change: not allowed for unit tracker, but possibly
            // OK for arrays on the boundary
            // TODO: disallow this and let higher-level TrackingGeometry
            // recover from the failure
            break;
        }

        state_ref_.cell    = init.cell;
        state_ref_.surface = init.surface;
        state_ref_.sense   = init.sense;

        // Find intersection
        Intersection next;
        try
        {
            next = tracker->intersect(this->state_ref());
        }
        catch (const std::exception& e)
        {
            ADD_FAILURE() << "Failed to find next surface at " << pos_
                          << " along " << dir_ << ": " << e.what() << endl;
            result.senses.push_back('x');
            result.distances.push_back(-inf);
            break;
        }

        result.distances.push_back(next.distance);

        // Move to surface
        if (next.distance == 0)
        {
            ++zero_movement;
            if (zero_movement >= max_zero_movements)
            {
                ADD_FAILURE() << "Encountered " << zero_movement
                              << " consecutive zero-distance movements at "
                              << pos_ << " along " << dir_ << endl;
                break;
            }
        }
        else if (next.distance != no_intersection())
        {
            zero_movement = 0;
            axpy(next.distance, dir_, pos_);
        }
        else
        {
            break;
        }
        state_ref_.surface = next.surface;
        state_ref_.sense   = next.sense;
    }
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Print result of a track
 */
void TrackerTest::TrackResult::print_expected() const
{
    cout << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
         << "static const char* const expected_cells[] = "
         << to_string(this->cells) << ";\n"
         << "EXPECT_VEC_EQ(expected_cells, result.cells);\n"
         << "static const char* const expected_surfaces[] = "
         << to_string(this->surfaces) << ";\n"
         << "EXPECT_VEC_EQ(expected_surfaces, result.surfaces);\n"
         << "static const int expected_senses[] = " << to_string(this->senses)
         << ";\n"
         << "EXPECT_VEC_EQ(expected_senses, result.senses);\n"
         << "static const real_type expected_distances[] = "
         << to_string(this->distances) << ";\n"
         << "EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);\n"
         << "/*** END CODE ***/\n";
}

//---------------------------------------------------------------------------//
} // namespace orange_test
