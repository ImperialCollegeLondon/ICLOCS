# Release Guide

This guide is for maintainers preparing a GitHub release.

## Pre-Release Checklist

- Confirm MATLAB CI passes on the release branch.
- Run representative examples locally with the intended solver stack.
- Regenerate `docs/API-INDEX.md` if MATLAB files changed:

  ```bash
  python3 tools/generate_api_index.py
  ```

- Review `README.md`, `docs/INSTALLATION.md`, and `docs/EXAMPLES.md`.
- Update `CITATION.cff` if version or citation details changed.
- Confirm the usage metrics workflow is healthy.

## GitHub Release Assets

Decision: future releases should attach explicit downloadable assets whenever maintainers want durable download counts.

GitHub only exposes download counts for release assets. Source-code ZIP downloads generated automatically by GitHub are not counted in the same way. If download tracking matters for a release, attach an explicit ZIP or MATLAB toolbox asset to the release.

Recommended assets:

- `ICLOCS-vX.Y.Z.zip`: a clean source distribution containing the repository files needed by MATLAB users.
- `ICLOCS-vX.Y.Z.mltbx`: optional MATLAB toolbox package if maintainers decide to support toolbox installation.

Suggested asset rules:

- Keep automatic GitHub source archives enabled, but do not rely on them for download metrics.
- Name release assets consistently so the usage metrics workflow can report cumulative release asset downloads.
- Do not include generated caches, solver binaries, local MATLAB preferences, or test output artifacts in release assets.
- Mention required external solvers and optional dependencies in the release notes rather than bundling third-party solver code.

## Suggested Release Notes

Include:

- Version number and date.
- MATLAB compatibility notes.
- Solver compatibility notes.
- Changes to examples, templates, or documentation.
- Known limitations.
- Migration notes for users updating existing problem files.
