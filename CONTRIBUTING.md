# Contributing To ICLOCS

Thank you for helping improve ICLOCS.

## Development Principles

- Preserve existing solver behavior unless a change explicitly targets solver behavior.
- Keep modernization work focused on tests, CI, documentation, examples, packaging, and repository metadata.
- Add tests before or alongside behavior changes.
- Avoid broad path additions that mix multiple example folders with duplicate function names.

## Before Opening A Pull Request

1. Run the MATLAB tests if MATLAB is available:

   ```matlab
   results = runtests('tests','IncludeSubfolders',true);
   assertSuccess(results)
   ```

2. Update documentation when user-facing behavior, setup steps, or examples change.
3. Keep pull requests focused and describe any solver/toolbox requirements.

## Reporting Issues

Please include:

- MATLAB version and operating system.
- Solver used, for example IPOPT, fmincon, or WORHP.
- Example or problem folder.
- Full error message and stack trace.
- Any local changes to `settings_*.m`.
