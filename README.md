# ORANGE (Celeritas-formatted derivative)

ORANGE (Oak Ridge Adaptable Nested Geometry Engine) is the new geometry engine
used for ray tracing and particle transport in
[SCALE](https://www.ornl.gov/scale). It has been declared non-export-controlled
as the visualization component of the [NEAMS
Workbench](https://www.ornl.gov/project/neams-workbench), as it has no
nuclear-specific components [0D999 classification](https://www.bis.doc.gov/index.php/documents/regulations-docs/2331-category-0-nuclear-materials-facilities-equipment-and-miscellaneous-items-1/file). The copyright on the code is still held by UT-Battelle (see LICENSE).

This repository is **not functional as a standalone code**. It serves as a
starting point for implementations in
[Celeritas](https://github.com/celeritas-project/celeritas) but generally is
not suitable for copying directly into the Celeritas codebase due to
infrastructure changes, low-level data refactoring, higher-level functional
refactoring, and stylistic inconsistencies.
