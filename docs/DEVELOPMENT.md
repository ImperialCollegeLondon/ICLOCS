# Development

The modernization work should preserve existing solver behavior. Prefer adding tests, documentation, CI, and repository metadata before changing MATLAB implementation files.

## Local Test Command

When MATLAB is available:

```matlab
results = runtests('tests','IncludeSubfolders',true);
assertSuccess(results)
```

The initial tests focus on pure helpers and problem-definition smoke checks. They intentionally avoid solving NLPs because solver availability varies across machines and GitHub runners.

## Continuous Integration

The MATLAB CI workflow runs on pull requests, pushes to `master` or `main`, and manual dispatch. It uses:

- `actions/checkout@v6`
- `matlab-actions/setup-matlab@v3`
- `matlab-actions/run-tests@v3`

The workflow writes JUnit test results and Cobertura coverage artifacts. Coverage currently focuses on files under `src`.

## Adding Tests

- Put MATLAB unit tests under `tests`.
- Keep fast unit tests separate from solver-heavy regression tests.
- Prefer tests that exercise public behavior rather than implementation details.
- Add solver-dependent tests later behind clear assumptions or separate workflows.

## Documentation Inventory

`docs/API-INDEX.md` inventories MATLAB files and their first help line. Use it to find files that still need richer help text.

## Website Migration

The original ICLOCS website contains substantial usage guidance. Track migration coverage in `docs/SOURCE_WEBSITE.md` and prefer moving durable setup, options, and example guidance into this repository so it can be versioned with the code.
