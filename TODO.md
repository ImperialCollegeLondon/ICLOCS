# ICLOCS Modernization TODO

This list tracks repository modernization work that should avoid changing the existing MATLAB solver behavior unless a future task explicitly calls for it.

## Priority 1: Usage Counter

- [x] Add README badges for recorded clone activity and release asset downloads.
- [x] Add a GitHub Actions workflow to refresh usage metrics.
- [x] Add a script that stores GitHub usage metrics in `docs/metrics/github-usage.json`.
- [x] Push the counter workflow and metric files to GitHub.
- [ ] Add a repository secret named `TRAFFIC_TOKEN` with repository Administration read access.
- [x] Run the **Update usage metrics** workflow manually once from the GitHub Actions tab.
- [x] Confirm the README badges update after the first successful workflow run.
- [x] Decide whether future releases should attach ZIP/toolbox assets so GitHub can count release asset downloads separately.

## Continuous Integration

- [x] Add a MATLAB-compatible CI workflow that runs on pull requests and pushes.
- [x] Check which MATLAB versions should be supported.
- [x] Add a lightweight syntax/path validation job that does not require expensive solver runs.
- [x] Add at least one smoke test that exercises a small example problem.
- [x] Add CI badges to the README once the workflow is stable.

## Unit Tests And Coverage

- [x] Identify stable public functions suitable for unit tests.
- [x] Add a MATLAB `runtests`-based test suite.
- [x] Start with low-cost tests around utilities, problem setup, interpolation, and input validation.
- [x] Add regression tests for a few known example problems without changing expected behavior.
- [x] Add code coverage reporting if the MATLAB CI environment supports it cleanly.

## Installation And Getting Started Documentation

- [x] Write first-time installation instructions for MATLAB users.
- [x] Document required MATLAB version, toolboxes, and external dependencies.
- [x] Add a quick-start example that gets a new user from clone/download to first successful solve.
- [x] Explain how to add ICLOCS to the MATLAB path.
- [x] Add troubleshooting notes for common setup and solver issues.

## Examples And Tutorials

- [x] Curate a small set of recommended beginner examples.
- [x] Add step-by-step walkthroughs for single-phase and multi-phase problems.
- [x] Document how to adapt a template into a new optimal-control problem.
- [x] Add expected outputs or screenshots/plots for key examples.
- [x] Mark advanced examples separately from beginner examples.

## Function And API Documentation

- [x] Inventory public functions, templates, and user-facing entry points.
- [x] Add or improve MATLAB help text for public functions.
- [x] Document inputs, outputs, defaults, side effects, and assumptions.
- [x] Generate browsable API documentation from MATLAB comments if practical.
- [ ] Keep internal implementation comments focused on behavior that is hard to infer from code.

## Repository Modernization

- [x] Add contribution guidelines.
- [x] Add issue and pull request templates.
- [x] Add a code of conduct if appropriate for the project.
- [x] Review repository layout and document what each top-level folder contains.
- [x] Add release guidance for maintainers.
- [x] Add citation metadata if there is a preferred paper or software citation.
