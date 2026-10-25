# OGS HM Tools

Educational notebook tools for custom OpenGeoSys HydroMechanics developments,
used alongside upstream OGSTools.

**Status: release preparation. No public installation package is available yet.**
The first prototype targets Linux x86_64 with CPython 3.12. It has been exercised
in a separate Python environment on the development host. A clean-machine test
and the complete bundled dependency license review remain pending.

## HM solver source

The custom OGS HM source is published on the
[`hm-development` branch](https://github.com/pavan0033/ogs-hm-tools/tree/hm-development/ogs).
See its [source guide](https://github.com/pavan0033/ogs-hm-tools/blob/hm-development/HM_SOURCE_README.md)
for contents, exclusions and build scope. The source snapshot is based on upstream
OGS plus the current HM changes; it includes custom unit tests but excludes large
upstream datasets, research results and the external imported HS material kernel.

Source commit for future toolkit integration:
[`f6fa80fd9ef1](https://github.com/pavan0033/ogs-hm-tools/commit/f6fa80fd9ef147e5ce4942687e01744e60c7b76d).
The prototype wheel predates this source publication and has not been rebuilt
from this commit. Public binary installation remains in preparation.

## Educational notebooks

- `notebooks/00_INSTALL_AND_CHECK_HM.ipynb`: install the supplied local prototype
  wheel, run a small MC HM example, and check computed pressure and hold fields.
- `notebooks/01_MCC_GEOSTATIC_WITH_SETUP.ipynb`: installation-enabled MCC teaching
  example. Upstream OGS runs its SMALL_DEFORMATION process and mesh utility;
  the custom HM runtime is installed alongside it for later HM examples.

Until a release is published, these notebooks require the prototype wheel from
the maintainer. Place it beside the notebook, or set `OGS_HM_WHEEL` to its path.
Their setup cells also install upstream Python dependencies from the internet.
Use a fresh Python 3.12 kernel. Installation requires Linux x86_64 for this build.

## Planned public installation

The reviewed wheel and its SHA-256 checksum will be attached to a versioned
GitHub Release. Notebook setup cells will then download that exact release.
This repository does not currently provide a working remote installation URL.
Users will not need access to the private development repository.

## Scope and validation

The HM prototype enables HydroMechanics; it does not replace upstream OGS for
other processes. The small column uses Constant-density incompressible fluid
and prescribed pressure at primary nodes. It is an installation demonstration,
not a boundary-driven seepage benchmark. Its 31 saved states passed pressure
and hold checks. A separate development audit checked weak mechanical balance,
MC admissibility and plastic history. The MCC notebook passed its original
checks for OCR 1, 2 and 5 in the same Python environment.

The imported HS kernel is excluded. Research results and commercial reference
images are excluded. See `release-status.json` for the exact candidate checksum
and outstanding release checks. This repository does not build the native solver;
installation will use prebuilt release assets when those are published.

## Credits and licensing

Based on OpenGeoSys, developed by the OpenGeoSys Community. The original OGS
license is retained in `licenses/OpenGeoSys.txt` and applies to OGS code.
The toolkit integrates with upstream OGSTools; it is not an official OGS or
OGSTools release. Contributor attribution for the HM additions remains to be
finalized by the author. No new blanket license for the author's contributions
or third-party bundled libraries is assigned by this preparation step.
