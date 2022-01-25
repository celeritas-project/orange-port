.. ############################################################################
.. File  : Geometria/orange/doc/orange.rst
.. ############################################################################

**************************************************
ORANGE: Oak Ridge Advanced Nuclear Geometry Engine
**************************************************

This chapter describes the implementation and key algorithms of ORANGE, a
constructive solid geometry (CSG) engine that implements particle tracking
with SCALE geometry definitions.

Introduction
============

The ORANGE implementation provides facilities for:

1. Construction: building a geometry model from processed SCALE input, from an
   ``.org.omn`` input format, or from other programmatic input.
2. Tracking: moving a particle through the geometry model, obtaining geometric
   cell, surface, and material IDs at each point.
3. Querying: converting the low-level geometry and tracking data (data
   structures and IDs) to meaningful user output.

The primary focus for the construction is extensibility, for tracking is
performance, and for querying is convenience.

Nomenclature
------------

Unit
   A KENO term for a combination of materials, other units, and arrays of other
   units, all enclosed by a boundary.

Universe
   A generalization of the unit that includes arrays. Particles can track along
   multiple universes simultaneously, with each universe corresponding to a
   particular depth or level. The units are related by transforms
   (rotations and translations) and inserted into one another via holes, array
   placement, and array elements.

Cell
   A volume of space with a unique identifier. Inside a universe, every point
   corresponds to exactly one cell. In a unit, cells are composed of a CSG
   tree of surfaces.

Surface
   Generally, surfaces are the boundaries between two cells. The same surface
   can be shared by multiple cell boundaries (e.g., the top part of a
   ``T``-shaped intersection).

Shape
   A shape is a utility for creating cells and surfaces. Each shape defines an
   *inside* region of space. This region is not necessarily convex.

Sense
   Whether a point or region is *inside* or *outside* relative to a
   surface or shape. Given a spherical surface with the formula
   :math:`x^2 + y^2 + z^2 - R^2 = 0`, the sense is *negative* when
   the expression on the left-hand side of the quadric equation is less than
   zero for a point and *positive* or when greater than or equal to zero. When
   referring to regions defined by shapes, *inside* and *negative* are
   equivalent, as are *outside* and *positive*. When evaluating a point in
   space, it can be logically *on* a surface when the quadric expression is
   approximately zero. Note that KENO uses ``+`` to signify *inside* (it
   includes or adds the shape) and ``-`` to signify *outside* (it excludes or
   removes the shape). ORANGE ``.org.omn`` input allows an MCNP-style *-/+*
   meaning
   *inside/outside* or a KENO-style boolean operation where omitting a sense
   means *inside* and prefixing with a tilde ``~`` means negating the shape for
   *outside*.

Region definition vector (RDV)
   A boolean expression given as the intersection of a list of pairs of senses
   and shapes. In KENO, an RDV defines each media, array placement, and unit
   boundary.

Proto
   The factory class for constructing units and arrays (shortened from
   *proto-universe*).

Level
   The depth below the global (top-level) unit of the universe that a
   particle is being tracked in. The following diagram explains the
   relationship between levels, universes, units, holes, and arrays.

State
   A complete description of how a certain particle relates to the geometry
   definition. In ORANGE, the state has *physical* components and *logical*
   components. The physical state is the global position and direction, which
   in theory define the particle's current cell and material. In practice, a
   logical state that comprises cell and surface information is needed to
   robustly track through the geometry.

.. figure:: levels.pdf
   :width: 500px

   Depiction of universes and levels. The circles are *units*, the square is an
   *array*, and the comments show how they are constructed from a KENO input.
   In this diagram, the same unit labeled *daughter* is used at two different
   levels.


Construction
============

Most of the complexity in construction is associated with the *unit proto*,
which is the equivalent of a KENO unit.

Deduplication
-------------

To improve the performance and robustness of tracking, it is often necessary to
simplify user input. Defining a geometry by closed 3D shapes typically means
that the constructing surfaces (planes, cylinders) are duplicated across
multiple shapes. Since the tracking within a unit in ORANGE is surface based
rather than volume based, it is necessary to combine equivalent (exactly or
nearly) surfaces into a single surface, so that crossing a surface means
leaving one shape and entering the next. This deduplication or elision of
"nearly equal" surfaces is controlled by a user-specified "fuzziness"
tolerance. This tolerance defaults to :math:`10^{-6}` for KENO geometry and
:math:`10^{-8}` for other geometry definitions.

Surface deduplication
^^^^^^^^^^^^^^^^^^^^^

When constructing surfaces, ORANGE applies local coordinate transforms (from
being inside a rotated shape) and then simplifies surfaces and their senses
into a canonical form. For example, a cylinder be defined parallel to the Z
axis and rotated to be parallel to the Y axis will first be transformed to a
*general quadric* surface (with ten unknowns), then simplified to a *simple
quadric* (with seven unknowns), then finally simplified to a Y-axis-aligned
cylinder (with three unknowns).

The region defined by a surface plus a sense (where the surface might be a
plane normal to the Z axis and the sense might be *negative* indicating
*below*) is also canonicalized during construction by flipping the sense and
the signs of the quadric surface coefficients. Currently this negation is only
applied to planes by ensuring that the first non-zero component along the x, y,
or z axis is in the positive half-space.

The final step of surface deduplication is to combine equivalent or
nearly-equivalent surface definitions to a single surface that can be reused
among cell definitions.  ORANGE does this by defining a hash-based set of
surfaces within a unit during construction, using a *soft* hash that truncates
floating-point values to a user-specified precision.

Shape deduplication
^^^^^^^^^^^^^^^^^^^

Additionally, users can define the same equivalent shape multiple times (and in
multiple different ways). Shape deduplication (more accurately, CSG node
deduplication) combines equivalent definitions for the same subregion of space.
This is most necessary for unit boundaries, which (with the exception of the
global unit) implicitly truncate daughter universes and are automatically
removed from daughter cell definitions.

During construction of a unit, ORANGE builds a CSG tree representation of the
input geometry. For all but the global unit, the boundary is locally implicit:
being in the unit means being inside its boundary, because boundary crossing
is controlled by the parent universe. This restriction implies restrictions on
other nodes of the CSG tree, which are simplified and eliminated from the rest
of the tree. For example, setting the boundary to a cuboid implies that every
interior object is to the right of the exterior ``-x`` surface, so that surface
is eliminated from every node in the tree. If the geometry description in the
output seems to have missing surfaces, it's probably because of this
deduplication.

Unit Proto
----------

Units are constructed with three different types of regions plus an external
boundary definition. The first region is a *media* entry, a cell filled with a
material. The second is an *array*, which is simply a daughter universe without
an explicit boundary, inserted into a local region (defined by local shapes).
The third is a *hole*, which is a daughter universe placed into the local unit.
Each hole entry has a *z-order* that defines whether the hole overrides other
regions in the unit. A z-order of *hole* means that the placed hole overrides
local media and arrays, but a *media* z-order can be specified to create holes
that are explicitly connected to other cells in the unit for a more
computationally efficient geometry.

Arrays and holes can have transforms, which are specified as the translation
needed to place the origin of the daughter universe to a particular point in
the local universe, a daughter-to-parent transform.

A unit's boundary is defined as the complement of a region definition
vector.  An option to make a KENO-like *implicit boundary* will cause the unit's
boundary to override local media, arrays, and holes. Otherwise, the boundary's
component shapes must be explicitly excluded from local media and holes.

.. note:: KENO restricts the global boundary RDV to shapes with "inside" senses
  only -- the intersections of the interiors of shapes -- but ORANGE lifts
  this restriction. Due to non-convex shapes it is possible to have reentrant
  boundaries even without the KENO restriction, so it is up to the analyst to
  ensure boundary conditions are physical and meaningful.

Array Proto
-----------

Array protos are 3-dimensional logically rectangular lattices of contiguous
units.
Three array types are implemented in ORANGE: rectangular, hexagonal, and
rhombic dodecahedral. Rectangular arrays support an arbitrary grid of cuboid
daughter cells, hexagonal arrays are triangular-pitched sets of hexagons, and
the dodecahedral arrays are the unit cell for face-centered-cubic lattices.

Arrays are defined with a 3-D array of daughter units. Each array type may
implement an "interior" function which allows hole-like placement in a parent
unit: that is, the array will create boundary shapes in the parent rather than
use existing boundary shapes. This capability is currently only implemented for
rectangular arrays, which can be placed as a cuboid in the parent geometry.

Hexagonal array protos support two additional input options that the others do
not. An input parameter can explicitly specify whether the array elements are
pointy-top or flat-top hexagonal prisms, though the default is to automatically
detect based on the input unit cells. A second input parameter allows the user
to input a *rectangular* grid of staggered hexagons for a more compact input
description.

Tracking
========

To enable the reuse of simple geometric constructs in multiple locations,
ORANGE tracks particles through multiple universes simultaneously. For
example, the top level might be the description of a reactor core, the next
might be an assembly, and the deepest level might be a pin cell. Each level can
have an independent coordinate system through transformations (translation and
rotation) applied when entering and leaving a level.

The ``Universe`` class is a thin composition of

 - A ``Tracker``, which stores geometric information such as surfaces and
   cells, and
 - A mapping of cells to daughter universes, which includes both pointers to
   the ``Universe`` and the corresponding transformations.

Each universe has a consecutive local set of unique *surface IDs* and *cell
IDs* representing distinct surfaces and cells in that level. Each cell is
comprised of bounding surfaces, which are *faces* of that cell.

A tracker's purpose is to find boundary crossings and logical (cell ID and
surface ID) state information for a particle through two primary methods:

- Initialize, which happens when sourcing a particle *or* when moving across a
  boundary; and
- Intersect, which calculates the next surface crossing (both distance and
  surface ID).

Multi-level tracking
--------------------

The ``Tracking_Geometry`` class is responsible for coordinating tracking across
the different levels of universes.  It also implements boundary conditions and
reflection, which are only applicable at the global level, and straight-lie
movement within a cell, which does not alter the logical state and therefore
does not require interaction with the trackers.

Initialization starts at the top-level universe (i.e., the global unit) and
asks the local tracker to initialize, yielding the local cell ID. If the cell
ID at the initialization level corresponds to a daughter universe, the geometry
state's level is incremented and the process repeats. When a non-universe cell
is encountered (i.e. a cell that corresponds to a physical material), the
process terminates.

Intersection (finding the closest boundary across all levels of tracking) is a
closed loop over the levels being tracked. The tracker at each level determines
the distance to the closest boundary. If the distance is *less than*
the nearest distance parent levels, that level is chosen as "closest". The
"soft" comparison (slightly expanding the apparent parent distance) effectively
makes a deeper-level distance look "a little more infinite" so that its
boundary is less preferable to that of the parent level. This
behavior is necessary to ensure that daughter universes only track in
their "interior" and never escape into an invalid region.

Intersection can only determine the next level to be "at or shallower" than the
current level: a deeper level can only be found after moving to an adjacent
cell at the lowest level.

Within-cell movement updates the distance to the next surface, but to improve
performance and reduce round-off error it does *not* update the local position
at any level. Positions are updated when moving the particle across a surface or
changing its direction.

Surface crossing is effectively the same as initialization, except that instead
of starting with a global coordinate, the deepest local coordinates (and
logical states such as exiting cell, crossing surface, and post-crossing sense)
are used to initialize.

Edge cases and robust tracking
------------------------------

Computational geometry for particle physics is notorious for tracking errors
caused by numerical imprecision. Errors typically occur at the boundary between
two regions. One solution to the imprecision is an exact bookkeeping of the
logical state of a particle by tracking its senses with respect to each
surface. This solution tends to fail when surface deduplication, multiple
universes, and cells with internal surface crossings are present. The obsolete
"GG" geometry implementation in SCALE attempted such bookkeeping at the cost of
expanding almost every instance of every unit and array unit in the geometry,
which was cost-prohibitive for many realistic geometries.

The basic requirement for avoiding edge cases is storing a logical state that
contains the current cell and surface. ORANGE stores cell and (TODO add this to
transport!) surface information at every level.

A robust geometry implementation is one in which no failures or infinite
loops are experienced for a well-defined geometry.
- Every point corresponds to exactly one region including a global exterior
  region.
- Each point can track to any other interior point in a finite number of
  surface crossings.

Additionally, ORANGE prohibits zero-distance steps as might occur
at edges or tangent surfaces between universes.

Edge cases
^^^^^^^^^^

Three key "edge" cases are present when tracking a particle along a straight
line:
- Initialization on a surface, either from a particle source (or raytrace
  origin) or from a higher-level universe. This initialization may occur
  exactly along a surface (and often does in ray-tracing and
  domain-decomposed applications).
- Finding the distance to the next cell boundary while on a surface. Such cases
  include those where the next surface *is* the current surface, e.g., moving
  across the inside of a cylinder.
- Crossing a surface, where the cell changes without a physical movement of
  the particle.

Another edge case appears during transport:
- A change in direction on or near a surface crossing. Direction changes
  *exactly* on a surface can occur in Shift when a user-defined boundary mesh
  specifies a reflecting boundary coincident with a surface in the problem.
  Direction changes very close to a surface may happen at the boundary between
  optically thin and thick regions of a problem.

An additional complication is that a particle may be physically "on" an
arbitrary number of surfaces: one when crossing an face; two when crossing an
edge or tangent surfaces; three or more when crossing vertices or tangent
surfaces.

One key feature of multi-level tracking in ORANGE is that the boundaries of
daughter universes are imposed by the parent universe, and daughters ideally
appear infinite in extent to avoid tracking into an invalid region (the
"exterior" of a daughter universe is not a valid region). The outer
boundaries of daughter units are automatically removed in most cases, which
eliminates the most common source of coincident boundaries (initialization on a
surface). However, unusual situations (including array placement inside a
shape that has surfaces coincident with that of the array daughter element's
shapes) can still cause multi-level coincident surfaces inside ORANGE.

Strategies
^^^^^^^^^^

Edge cases on universe boundaries need to be considered in multi-level
tracking.  This class of edge case is about finding the "correct"
universe in which to initialize a particle, while accounting for the
possibility of tiny regions of "invalid" space on and around the outer boundary
of the daughter universe. A non-robust implementation could result in an
infinite loop or a geometry error (no valid region) when crossing or
initializing on problematic boundaries.  There are two broad ways to handle
boundary error mitigation.

The first approach guarantees that universe
locating will converge: when descending into successive daughter universes,
each daughter *must* successfully initialize the particle if it is within a
certain tolerance of the valid region. Since initialization cannot change the
physical position of particles, the logical state may not exactly agree with
the physical state (e.g. a particle may be logically "inside" one or more
surfaces needed to place it inside a valid region, while still having a
physical position on the wrong side of that surface). The potential
inconsistency in the state may lead to edge cases in distance-to-boundary
evaluations and surface crossings. Although this approach can guarantee
initialization, searching for a valid *inside* logical state given a slightly
out-of-bounds position requires heuristics that are themselves error-prone. For
example, if initialization fails to find a valid logical state, one can try
evaluating a nearby physical position (temporarily "bumping" the particle), and
there are numerous ways to guess at nearby valid locations, none of which is
bulletproof.

.. note:: As of this development version, (TODO update this later) the current
   approach in ORANGE is the first approach. The array trackers have extra
   logic to try to force particles on or just outside the edge of the array
   boundaries into a valid array unit. The unit trackers have even more
   complicated logic. If the particle is "on" (crossing) a surface, the logic
   evaluator attempts both senses with respect to the current face. The masked
   tracker has to prioritize whether the sense was flipped or not in order to
   prevent a surface from appearing to be inside multiple cells.

   On one branch of ORANGE I tried testing for being exactly on a surface, as
   this was a common failure case for raytrace initialization. This would allow
   initialization on a surface even from the start of a ray.

   None of the approaches (or even the more rigorous GG approach) solve the
   case of initialization on points that are exactly on the boundary of
   multiple surfaces (such as the corner of a hex).

A second approach is to allow universes to reject invalid regions and bump the
particle as part of the universe location loop if
no universe can correctly initialize the particle at a position. This allows
universes to be very strict about which positions are valid, ensuring
consistency between logical and physical states. The downside is the potential
of infinite or arbitrarily-long loops (e.g., bumping the particle
:math:`O(1/\epsilon)` times given a geometry tolerance `\epsilon`).

Invalid regions may include:

  * A coincident surface between a higher and lower level that was not
    eliminated through surface deduplication. (Since the position of the
    surface on multiple levels might be slightly different due to floating
    point errors in transformation *and* might be shifted due to surface
    deduplication, the region enclosed by this invalid "surface" may actually
    have a width of :math:`O(\epsilon)` along its thinnest dimension rather
    than being a true surface.)
  * The exact edges or corners where two or more surfaces intersect.

The updated logic in ORANGE is:

- Try to initialize (from global position or from surface crossing) at the
  original position. Let universes reject invalid positions.
- If no universe accepted the point as valid, bump first along the direction of
  movement and try initializing again.
- If initialization still failed, try at the corners of an axis-aligned
  cube (with a circumradius of the bump distance) that do not move the particle
  backward. (In other words, the dot product of the bump vector and the
  direction must be nonnegative.)

If the bump results in a component orthogonal to the particle's direction of
movement (i.e. the particle is traveling exactly along a cartesian face) then
that component is added to a cumulative "correction" term that is carried
around for the lifetime of the particle. The correction term is only used to
report the particle's position back to the rest of the Shift code base and is
not used for internal tracking. It will be reset if the state is reinitialized.
Future work could use this component and a "safety" distance (distance to the
nearest boundary in *any* direction) to allow efficient movement within the
cell in any direction, a capability needed for transport of charged particles
in electric or magnetic fields.

In order to allow the calculation of surface normals while tracking, ORANGE
needs to know the ID of a surface, its universe, and the current cell for every
universe. This latter requirement is because each cell can have a different
transformation to the parent, and surface normal calculations are always
performed in the lowest frame of reference, and ORANGE does not cache
local-to-global transformations.

Bumping
^^^^^^^

The distance of a "bump" in ORANGE is determined by the user-specified global
"fuzziness" tolerance, the characteristic length scale (which defaults to 1
cm), and the
particle's local position. The bump distance must be large enough to guarantee
crossing surfaces that were effectively "moved" during surface deduplication
(where the location of the surface in the parent universe is slightly different
than the location in the daughter universe), and large enough to ensure the
movement is greater than the limits imposed by floating point precision (e.g.,
adding a fixed value of :math:`10^{-12}` to a coordinate with a value of 1000
will not change the floating-point representation in double-precision
arithmetic).

Different heuristics useful for bumping (both physically and logically)
include testing one or more positions:
- along the direction of particle's travel,
- on the faces or corners of a small cube, either with axis-aligned faces or
  with a face normal along the direction of the particle's travel,
- toward the center of the local universe, and
- along the inward normal toward the closest surface.

Simple Unit Tracker
-------------------

The simple unit tracker is a fast implementation for units that do not have
any overlapping regions: all components are directly connected to surfaces. A
cell is still allowed to truncate a descendant universe; it simply cannot
overlap with other cells in the same unit. This is the unit type created by
KENO-V.a geometry, VERA geometry, Oberon, etc. It does *not* support
overlapping cell detection.

Initialization
^^^^^^^^^^^^^^

Initialization is implemented by searching for the spatial position in a
bounding volume hierarchy (BVH) or other acceleration structure to find a
subset of the unit's regions that might contain the point. For each candidate
cell, the routine evaluates the surface senses of the cell's faces, then checks
whether the cell's logic places the particle inside. The first "inside" result
stops the search for the sake of efficiency, although this means that "simple
units" do not check for overlapping volumes.

When crossing a cell, the starting cell is skipped in the loop over candidate
cells. If the check of surfaces indicates the point is exactly coincident with
a surface but was not crossing that surface, initialization returns a failure
so that the multi-level tracking logic will bump the particle off the surface
in the interest of robustness.

Intersection
^^^^^^^^^^^^

There are two algorithms for finding the nearest surface crossing. The first is
the "simple" one, where crossing any boundary exits the cell. In that case,
the closest nonzero distance points to the target surface. The second "complex"
algorithm is for cells with internal surfaces. The determination of whether a
cell is "complex" is done during construction and is currently can be
conservative, e.g. by declaring that only regions comprised of intersections
between convex shapes are considered "simple." A more accurate but much
more expensive determination is to use the CSG definition of a cell can be
evaluated over all possible combination of surface states: if any sense flip
when the cell is in an "inside" state does *not* result in a change of
inside/outside, then the cell is complex. (This may still be unnecessarily
conservative if some of the combination of senses is not physically possible.)

The "complex" algorithm effectively tracks through the internal surfaces in
order until outside the cell:
- Calculate all senses at the current point, accounting for a potential current
  surface ID and sense used when finding the distance immediately after a
  boundary crossing
- Calculate intersections for all surfaces in the cell
- Partition the found intersections to separate strictly positive finite
  distances from invalid or infinite distances
- Sort the potential intersections in ascending order by distance-to-intercept
- Loop over each surface crossing:
  * Flip the calculated sense for the face being crossed
  * Evaluate the cell logic expression with the updated senses
  * If the evaluation is now "outside" the cell, then the next surface and
    sense are saved and the search is complete.

Acceleration structures
^^^^^^^^^^^^^^^^^^^^^^^

The methods to accelerate initialization and intersection rely on reducing the
search space over cells and surfaces from linear to logarithmic or better.
Under the right conditions, the cells within a unit can be located in
logarithmic time by constructing a bounding hierarchy such as a kd-tree or
octree, or constant time by constructing a fixed-width mesh. In this section,
any of these acceleration structures will be referred to as a "grid" and the
smallest element of the structure a "grid cell."

Accelerating initialization is easy.  Each grid cell of the acceleration
structure must contain a list of geometry cells that may intersect it. This
list of cells is determined with the help of axis-aligned bounding boxes
constructed by each of the shapes used to specify the cell.  The acceleration
algorithm locates the grid cell containing the given point and thus obtains a
list of all possible cells to loop over.

Accelerating *intersection* on the other hand is not easy. It requires
initializing the grid coordinates for the particle in a local state at the
start of the intersection test, then tracking along the acceleration grid until
the surface crossing is found.

The "simple" cell intersection method progressively tests all potential
surfaces in each grid cell. When an intersection is found, that surface
crossing is the result. Determining the surfaces that pass through a grid cell
is problematic in the general case of arbitrary quadrics. In the case of the
simple intersection algorithm, the non-conservative approach is to assume that
the surfaces of any volume may be in any grid cell containing the volume. This
is acceptable because in the simple method, a surface crossing happens if and
only if the volume changes.

The accelerated "complex" cell intersection algorithm is similar to the
non-accelerated case,
but instead of evaluating intersections for all grid cells and considering
distances :math:`0 < d < \infty`, each grid cell evaluates only surfaces that
that intersect the grid cell and then flips the bits only for those whose
intersection distances are :math:`d_\mathrm{i} < d \le
d_\mathrm{f}` where the bounds are those of the grid cell's intersection with
the particle's track. Unlike the simple case, it is necessary to actually
calculate intersections between the grid cells and the surfaces for the complex
tracking scheme.

The initialization scheme will never be much slower with acceleration, assuming
the cost of finding the starting point in the grid is cheap, because the list
of cells to check will always be smaller than the list of all cells in the
unit. However, intersection testing could conceivably be much slower for
complicated units or those with cells that enclose many grid cells, because the
same bounding surfaces might have to be tested multiple times.

Masked Unit Tracker
-------------------

The "masked unit tracker" is for tracking particles in a KENO-VI unit. At a
high level, a unit comprises:

 - Shapes, which themselves comprise one or more quadric surfaces;
 - Media, specified with a Region Definition Vector (RDV) of boolean operators
   on shapes;
 - Arrays, which are inserted into a unit by translating an array definition
   into the local unit and filling a cell of space specified by
   an RDV;
 - Holes, which are other units that "mask" (replace or take precedence over)
   media and arrays in this unit; and
 - A boundary, also defined with an RDV.

ORANGE internals translate these definitions into *surfaces*, which are
low-level quadric definitions, and *cells*, which specify potentially
overlapping regions in space. Each cell is a boolean expression of surfaces
that defines a region in space, and each cell has a "Z order" that specifies
precedence over other cells. The boundary of a unit (a cell defining the
exterior) has the highest Z order, followed by holes, then arrays, then media.
This allows user-specified regions to implicitly truncate others.

Initialization
^^^^^^^^^^^^^^

Initialization in a unit tracker uses a kd-tree acceleration structure to find
a list of all possible local cells that enclose the point. These cells are
iterated over in descending Z order. A cell is "found" by calculating the
senses of the point with respect to the cell's enclosing faces and evaluating
the cell's "inside" logic expression with those senses. When the particle is on
a surface, the sense with respect to the surface is known a priori, so
the logical state of a particle can change from one side of a surface to
another without any change in its physical position. If the particle is
crossing a surface, the cell it starts in is skipped, as particles cannot cross
a surface and end up in the starting cell.

To ensure that the particle's new location is not simultaneously in multiple
cells (not in overlapping regions), the search continues for other potentially
"inside" cells. When a lower Z order than the first "found" cell is
encountered, then the search is complete (since all less important cells are
being hidden by the higher Z order cells).

Intersection
^^^^^^^^^^^^

The intersection for a masked tracker has to be cell-based rather than
surface-based, because the z-order is more important than the surface
connectivity. However, we can use the surface/cell connectivity to build lists
of possibly important surface crossings at construction. ::

   LOOP over z-order less than or equal to exiting cell
      CALC distance to surfaces with this z-order
      LOOP over surfaces closer than current min position
         CALC updated particle position just past this surface
         LOOP cells connected to this surface with current z-order
            CALC surface senses for the cell at this surface crossing
            IF inside a new cell (or outside the starting cell)
               UPDATE the new minimum distance, surface, and sense
               BREAK out of this loop, check next surface


Rectangular Array Tracker
-------------------------

The rectangular array tracker is an implementation detail. It internally
indexes with "row-major" cell IDs: the ``x`` axis varies the fastest, and
``z`` the slowest. The surface IDs are the same as those used in Cartesian
nodal tallies (see the ``surface_indexing`` technical note).

To allow for nearly-coincident surfaces between arrays and the enclosing shapes
or unit, any point spatially outside the array tracker is still *logically*
inside the array. The array boundaries are "bumped" outward to account for
this.

Initialization
^^^^^^^^^^^^^^

Intersection
^^^^^^^^^^^^

Querying
========

The high-level metadata classes correspond to low-level trackers.

Shift interaction
=================

The top-level ``ORANGE_Geometry`` class handles the interaction of ORANGE state
information and the Shift geometry tracking API. Because Shift will make
repeated calls to geometry information that might not have changed (position,
next-distance, etc.), the geometry class calculates caches such data as
necessary. Cell IDs (which require local-to-global ID mapping) and materials
are calculated on demand. The call to cross surfaces, find
distance-to-boundary, etc. are dispatched to the trackers. Queries for cell
and surface names are dispatched to the Metadata classes.

Overlapping region errors caught in the Masked Unit Tracker by crossing
surfaces or initializing particles are caught by the Shift interface. It
constructs a message using the metadata for the failing universe, including
the cell names, the *local* (within-unit) cell IDs, and the line number where
the cell was defined. The local-to-global translation can
be performed by looking at the geometry's description, finding the failing
unit, and adding the local cell ID to that universe's cell offset.


.. ############################################################################
.. end of Geometria/orange/doc/orange.rst
.. ############################################################################
