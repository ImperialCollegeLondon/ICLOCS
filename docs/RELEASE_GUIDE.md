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

GitHub only exposes download counts for release assets. Source-code ZIP downloads generated automatically by GitHub are not counted in the same way. If download tracking matters for a release, attach an explicit ZIP or MATLAB toolbox asset to the release.

## Suggested Release Notes

Include:

- Version number and date.
- MATLAB compatibility notes.
- Solver compatibility notes.
- Changes to examples, templates, or documentation.
- Known limitations.
- Migration notes for users updating existing problem files.
