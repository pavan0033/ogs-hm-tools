# Custom OGS HydroMechanics source

The `hm-development` branch contains the custom OGS source under `ogs/`.
The `main` branch contains the educational toolkit landing page and notebooks.

This is a working-tree snapshot based on upstream OGS commit
`c02fe15b294d3acd576a8d6e1736f0a38ab06cd8`, including the uncommitted HM additions.
It imports no private development Git history. `HM_SOURCE_MANIFEST.json` records
all exported hashes and the lists of changed and added source files.

## Contents

The source includes HM engineering phases, initialisation, hydraulic controls,
excavation/construction, checkpoint/restart handling, safety continuation,
structural coupling, transactional material integration and solver controls.
It also includes the custom material/host unit tests. Existing upstream copyright
and license notices are retained. OGS licensing is in `ogs/LICENSE.txt`.

The external imported PublishedHS material kernel is not included. Generic
integration hooks and tests referring to that optional behaviour remain in the
source; tests requiring its external library skip when it is not supplied.
The upstream numerical benchmark datasets (`Tests/Data`) and website/tutorial content
(`web/`, except build metadata in `web/data/`) are omitted to keep this source branch focused. They can be
obtained separately from the upstream revision when needed. Local research
outputs, commercial reference images and runtime binaries are not included.

## Build and validation scope

Follow the compiler, CMake, Ninja and dependency prerequisites in `ogs/README.md`
and `ogs/CMakePresets.json`. This source snapshot requires C++23 and CMake 3.31+
according to the current source instructions. MFront-enabled builds require the
matching TFEL/MGIS dependencies. Do not assume the toolkit's binary wheel or
upstream pip package reproduces these source changes automatically.

For an initial solver-only build, use an isolated build directory, enable the
HydroMechanics process and MFront, and disable OGS_BUILD_TESTING and OGS_BUILD_UTILS.
The omitted upstream test data are needed for the full upstream benchmark suite.
The source is published without rebuilding or modifying the active development
runtime. Export validation verifies exact file hashes; it is not a new compile,
a new numerical qualification, or a clean-machine portability result.

## Future toolkit link

A toolkit release should record the exact commit of this branch plus build
configuration and dependency versions, then use binaries built from that commit.
Pin a full commit or a versioned tag rather than a moving branch in released
notebooks. Do not represent this as an official upstream OGS release.

Binary distribution license review and clean-machine tests remain separate from
this source publication. Current qualified feature scopes and pending studies
must be checked before making numerical capability claims.
