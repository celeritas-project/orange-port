//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/DodeArrayTracker.hh
 * \brief DodeArrayTracker class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Tracker.hh"

#include "base/RegularIndexer.hh"
#include "orange/surfaces/SurfaceContainer.hh"

namespace nemesis
{
template<class T>
struct RhombicDodecahedronFaceTraits;
}

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Rhombic dodecahedral array tracker.
 *
 * The rhombic dodecahedral tracker supports pointy-topped (+Z) rhombic
 * dodecahedrons (hexagonal XZ/YZ cross section). The dodecahedrons are
 * arranged into a rectangular vertical and square horizontal layout.
 *
 * The dodecahedron is described by an inscribing radius, or apothem, making
 * the horizontal square pitch equal to 2*apothem and the vertical rectangular
 * height equal to sqrt(2)*2*apothem.
 *
 * Internally the 12 faces are described using x,y and w,i,j,k axes. The
 * w,i,j,k axes describe the 4 oblique planes that constitute the pointy-top of
 * the dodecahedron. The oblique planes present the following projection
 * from +Z looking into -Z.
 *
 * <pre>
 *   (+Y)                             (+Y)
 *    ^      [ Top ]                   ^      [ Bottom ]
 *    |                                |
 *    |   +--------+--------+          |   +--------+--------+
 *    |   |        |        |          |   |        |        |
 *    |   |   K    |   W    |          |   |  -I    |  -J    |
 *    |   |(-x+y+z)|(+x+y+z)|          |   |(-x+y-z)|(+x+y-z)|
 *    |   +--------+--------+          |   +--------+--------+
 *    |   |        |        |          |   |        |        |
 *    |   |   J    |   I    |          |   |  -W    |  -K    |
 *    |   |(-x-y+z)|(+x-y+z)|          |   |(-x-y-z)|(+x-y-z)|
 *    |   +--------+--------+          |   +--------+--------+
 *    +-----------------------> (+X)   +-----------------------> (+X)
 * </pre>
 *
 * The surface number, axes, and normals are described in the table below
 *
 *   | Surf | Axis | Normal                        |
 *   |------|------|-------------------------------|
 *   | 0    | -X   | <-1,0,0>                      |
 *   | 1    | +X   | <1,0,0>                       |
 *   | 2    | -Y   | <0,-1,0>                      |
 *   | 3    | +Y   | <0,1,0>                       |
 *   | 4    | -W   | <-half,-half,-sqrt(2)*apothem>  |
 *   | 5    | +W   | < half, half, sqrt(2)*apothem>  |
 *   | 6    | -K   | < half,-half,-sqrt(2)*apothem>  |
 *   | 7    | +K   | <-half, half, sqrt(2)*apothem>  |
 *   | 8    | -I   | <-half, half,-sqrt(2)*apothem>  |
 *   | 9    | +I   | < half,-half, sqrt(2)*apothem>  |
 *   | 10   | -J   | < half, half,-sqrt(2)*apothem>  |
 *   | 11   | +J   | <-half,-half, sqrt(2)*apothem>  |
 *
 *  Note: The axis vectors can be flipped via XOR bitwise operator.
 *        Also, this assumes the surface ordering specified in
 *        the RhombicDodecahedronShape
 *
 * XY cross section at even and odd Z indices are below. The odd z indices
 *  are displaced in XY by the length of the dodecahedron's apothem.
 *
 * <pre>
 * (+Y)                                (+Y)
 *  ^      [ Z = 0 ]                    ^      [ Z = 1 ]
 *  |                                   |
 *  |                                   |
 *  |                                   |       |--------|--------|
 *  |                                   |       |        |        |
 *  |   |--------|--------|             |       |0,1,1   |1,1,1   |
 *  |   |        |        |             |       |        |        |
 *  |   |0,1,0   |1,1,0   |             |       |--------|--------|
 *  |   |        |        |             |       |        |        |
 *  |   |--------|--------|             |       |0,0,1   |1,0,1   |
 *  |   |        |        |             |       |        |        |
 *  |   |0,0,0   |1,0,0   |             |       |--------|--------|
 *  |   |        |        |             |
 *  |   |--------|--------|             |
 *  +-----------------------> (+X)      +-----------------------> (+X)
 * </pre>
 *
 * XZ cross section depicting W,I,J,K, and X,Y axes index change as a function
 *  of the Z index.
 *
 * clang-format off
 *
 * <pre>
 * (+Z)                                     (+Z)
 *  ^      [ Y = 0 ]                         ^      [ Y = apothem ]
 *  |                                        |
 *  |                                        |   \ 0,0,4/ \ 1,1,4/ \ 2,1,4/ \
 *  |                                        |     \ /     \ /      \ /      \
 *  |                                        |      |        |        |       |
 *  | | |0,0,3   |1,0,3   |2,0,3   |         |       / \ 0,0,3/ \ 1,0,3/ \
 * 2,0,3/     | | | |        | |    /      \ /      \ /      \  /       |    /
 *     \     /  \     / \     / |   |        |        |        |         |  /
 * 0,0,2\ / 1,0,2\ / 2,0,2\ /
 *  |   |0,0,2   |1,0,2   |2,0,2   |         |  \       / \      / \      / \
 *  |   |        |        |        |         |    \  /K    W\ /      \ /      \
 *  |    \     /  \-J -K/  \     /           |      |        |        | | | \
 * /0,0,1\ / 1,0,1\ / 2,0,1       |      |0,0,1   |1,0,1   |2,0,1   | | / \ /
 *            \      / \            |      |        |        |        | |    /K
 * W\ / \ / \          |    /  \     /  \     /  \     / |   | | | | |  /0,0,0
 * \ /
 * 1,0,0\ / 2,0,0\ / |   |0,0,0   |1,0,0   |2,0,0   |         | \     /  \ /  \
 * / |   |        |        |        |         |    \ / \ /      \ / |    \ /  \
 * /  \     /           | |      \ /      \ / \ /             |
 *  |--------------------------------> (+X)  |---------------------------->
 * (+X)
 * </pre>
 *
 * clang-format on
 *
 *  The index scheme is as follows:
 *
 * | Surf | Axis | U,V,W (Even Z) | U,V,W (Odd Z) | Index Transform Vector |
 * |------|------|----------------|---------------|------------------------|
 * | 0    | -X   | <-1, 0, 0>     | <-1, 0, 0>    | null                   |
 * | 1    | +X   | < 1, 0, 0>     | < 1, 0, 0>    | null                   |
 * | 2    | -Y   | < 0,-1, 0>     | < 0,-1, 0>    | null                   |
 * | 3    | +Y   | < 0, 1, 0>     | < 0, 1, 0>    | null                   |
 * | 4    | -W   | < 0, 0,-1>     | <-1,-1,-1>    | -< 1, 1, 0>            |
 * | 5    | +W   | < 0, 0, 1>     | < 1, 1, 1>    |  < 1, 1, 0>            |
 * | 6    | -K   | < 1, 0,-1>     | < 0,-1,-1>    | -< 1, 1, 0>            |
 * | 7    | +K   | <-1, 0, 1>     | < 0, 1, 1>    |  < 1, 1, 0>            |
 * | 8    | -I   | < 0, 1, -1>    | <-1, 0,-1>    | -< 1, 1, 0>            |
 * | 9    | +I   | < 0,-1, 1>     | < 1, 0, 1>    |  < 1, 1, 0>            |
 * | 10   | -J   | < 1, 1,-1>     | < 0, 0,-1>    | -< 1, 1, 0>            |
 * | 11   | +J   | <-1,-1, 1>     | < 0, 0, 1>    |  < 1, 1, 0>            |
 *
 * Note: the index vector U,V,W changes as function of Z being even or odd.
 *       This is because the rhombic dodecahedral array is rectangular
 *       vertically instead of triangular. However, the Index transform vector
 *       < 1, 1, 0> accounts for this.
 */
class DodeArrayTracker final : public Tracker
{
  public:
    //@{
    //! Public type aliases
    using DimVector   = Array<size_type, 3>;
    using SenseVector = Array<Sense, 12>;
    using Face        = ArrayFace<RhombicDodecahedronFaceTraits>;
    //@}

  public:
    //// CONSTRUCTION ////

    // Construct from rhombic dodecahedral mesh
    explicit DodeArrayTracker(real_type apothem, const DimVector& dims);

    //// TRACKING ////

    // Find the local cell and possibly surface ID.
    Initialization initialize(LocalState state) const final;

    // Calculate distance-to-intercept for the next surface
    Intersection intersect(LocalState state) const final;

    // Calculate normal on the current surface
    Real3 normal(LocalState state) const final;

    //// ACCESSORS ////

    //! Number of cells
    size_type num_volumes() const final { return cell_indexer_.size(); }

    //! Number of surfaces
    size_type num_surfaces() const final { return surfaces_.size(); }

    //! Dimensions
    const DimVector& dims() const { return cell_indexer_.dims(); }

    //! Cell id for given coordinate
    VolumeId volume_id(const DimVector& c) const { return coord_to_cell(c); }

    //! Dodecahedron apothem (interior radius)
    real_type apothem() const { return apothem_; }

    // Find the cell coordinates containing the point
    DimVector find(SpanConstReal3 point) const;

    // Access all surfaces
    const SurfaceContainer& surfaces() const { return surfaces_; }

  private:
    //// TYPES ////

    using CellIndexer_t = RegularIndexer<size_type, 3>;

    //// DATA ////

    // Stores indices as [X, Y, Z]
    CellIndexer_t cell_indexer_;

    // Interior radius
    real_type apothem_;

    // Vector of surface definitions for rhombic dodecahedron
    SurfaceContainer surfaces_;

    // Vector of "inside" surface sense for each face
    SenseVector senses_;

    //// IMPLEMENTATION ////

    DimVector cell_to_coords(VolumeId cell) const;
    VolumeId  coord_to_cell(const DimVector& ijk) const;

    Initialization initialize_interior(LocalState) const;
    Initialization cross_surface(DimVector, Face) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
