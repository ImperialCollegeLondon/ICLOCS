# Public API Reference

This page documents the MATLAB entry points most users are expected to call directly. For a file-by-file inventory of solver internals, generated helper functions, and example callbacks, see [`API-INDEX.md`](API-INDEX.md).

ICLOCS problem definitions follow the same pattern:

```matlab
[problem, guess] = MyProblem;
options = problem.settings(100);
[solution, status] = solveMyProblem(problem, guess, options);
```

## Problem Definition Functions

Each example problem has a top-level function named after the example, such as `BangBang`, `DoubleIntegratorTracking`, or `OrbitRaising`.

Syntax:

```matlab
[problem, guess] = ProblemName;
```

Inputs:

- None.

Outputs:

- `problem`: structure containing time bounds, state bounds, input bounds, path and boundary constraint bounds, solver callbacks, scaling data, and user data.
- `guess`: structure containing initial trajectories for states, inputs, parameters, and time.

Side effects and assumptions:

- These functions do not solve the NLP. They only define the optimal-control problem and initial guess.
- The problem folder must be on the MATLAB path so related settings, dynamics, and derivative callbacks can be found.
- The returned `problem.settings` field is a function handle to the example-specific settings file.

## Settings Functions

Each problem normally provides a `settings_*` function.

Syntax:

```matlab
options = settings_ProblemName(n);
options = settings_ProblemName;
```

Inputs:

- `n`: optional mesh-density or discretization-size value used by many examples. The exact meaning is example-specific.

Outputs:

- `options`: solver and transcription options, including discretization, mesh strategy, solver selection, tolerances, scaling, derivative settings, and print/plot controls.

Defaults:

- When `n` is omitted, each example chooses its own default mesh or discretization size.
- Most examples are configured for direct collocation with IPOPT when IPOPT is available.

Side effects and assumptions:

- Settings functions construct an options structure only. They do not run the solver.
- Solver-specific options require the corresponding solver to be installed and on the MATLAB path.

## `solveMyProblem`

Main entry point for transcribing and solving an ICLOCS problem.

Syntax:

```matlab
[solution, status] = solveMyProblem(problem, guess, options);
[solution, status, OCP] = solveMyProblem(problem, guess, options);
[solution, history] = solveMyProblem(problem, guess, options);
[solution, history, OCP_MR, OCP_initial] = solveMyProblem(problem, guess, options);
[solution, status] = solveMyProblem(problem, guess, options, OCP_in);
```

Inputs:

- `problem`: problem structure returned by a problem-definition function.
- `guess`: initial-guess structure returned by a problem-definition function or a warm-start routine.
- `options`: options structure returned by the problem settings function.
- `OCP_in`: optional pre-transcribed problem structure from an earlier call, used when re-solving an updated problem without repeating transcription.

Outputs:

- `solution`: post-processed solution structure. Typical fields include time, state, input, cost, interpolation polynomials, local error estimates, constraint violations, and computation time.
- `status`: solver status for fixed-mesh solves.
- `history`: mesh-refinement history for mesh-refinement solves, including solution, status, error, constraint-error, and timing histories.
- `OCP`, `OCP_MR`, `OCP_initial`: transcribed NLP data structures that can be reused for warm starts or problem updates.

Defaults:

- Defaults are supplied through `options`, not inside `solveMyProblem`.
- Post-solve analysis runs unless the relevant `skipPostAnalysis` option is enabled.

Side effects and assumptions:

- Clears intermediate ICLOCS variables before solving.
- Calls the selected NLP solver and may create figures or print progress depending on `options.print`, `options.plot`, and multiphase equivalents.
- Requires solver dependencies, problem callbacks, and ICLOCS `src` folders on the MATLAB path.
- Mesh-refinement mode may run multiple NLP solves and return refinement history instead of a simple solver-status structure.

## `updateMyProblem`

Update selected fields of a previously transcribed single-phase direct-collocation problem.

Syntax:

```matlab
OCP_updated = updateMyProblem(OCP_old, 'x0', x0);
OCP_updated = updateMyProblem(OCP_old, 'z0', z0);
OCP_updated = updateMyProblem(OCP_old, 'bl', bl, 'bu', bu);
OCP_updated = updateMyProblem(OCP_old, 'userdata', userdata);
```

Inputs:

- `OCP_old`: transcribed OCP structure returned by `solveMyProblem` when the appropriate third output is requested.
- Name-value updates:
  - `'x0'`: new initial state vector.
  - `'z0'`: new NLP initial decision vector.
  - `'bl'`: new boundary-constraint lower bounds.
  - `'bu'`: new boundary-constraint upper bounds.
  - `'userdata'`: structure merged into the existing user data.

Outputs:

- `OCP_updated`: updated transcribed OCP structure suitable for a subsequent `solveMyProblem(..., OCP_updated)` call.

Defaults:

- Omitted name-value pairs are left unchanged.

Side effects and assumptions:

- Supports single-phase direct-collocation problems only.
- Applies scaling to `x0` when scaling is enabled in the stored OCP data.
- Does not solve the updated problem by itself.

## `simulateSolution`

Simulate an optimized open-loop solution with MATLAB ODE integration.

Syntax:

```matlab
[tv, xv, uv] = simulateSolution(problem, solution);
[tv, xv, uv] = simulateSolution(problem, solution, odesolver);
[tv, xv, uv] = simulateSolution(problem, solution, odesolver, tstep);
[tv, xv, uv] = simulateSolution(problem, solution, odesolver, tstep, stateInputSwap);
```

Inputs:

- `problem`: original problem structure.
- `solution`: solution returned by `solveMyProblem`.
- `odesolver`: optional solver name, one of `'ode113'`, `'ode45'`, or `'ode23'`. The default is `'ode113'`.
- `tstep`: optional scalar time step used to build `solution.t0:tstep:solution.tf`.
- `stateInputSwap`: optional two-row index map used by examples that swap selected states and inputs for simulation.

Outputs:

- `tv`: simulation time vector returned by the ODE solver.
- `xv`: simulated state values.
- `uv`: interpolated input values evaluated at `tv`.

Side effects and assumptions:

- Uses `problem.sim` callbacks and interpolation data stored in `solution`.
- Clips simulated/interpolated inputs to the configured input and mapped-state bounds.
- Simulates the supplied solution; it does not re-optimize.

## `simulateSolutionSegment`

Segment-oriented variant of `simulateSolution` for workflows that provide an explicit time vector.

Syntax:

```matlab
[tv, xv, uv] = simulateSolutionSegment(problem, solution, odesolver, T);
[tv, xv, uv] = simulateSolutionSegment(problem, solution, odesolver, T, stateInputSwap);
```

Inputs and outputs:

- Same as `simulateSolution`, except `T` is the complete integration time vector rather than a scalar time step.

Side effects and assumptions:

- Uses the same supported ODE solvers and interpolation behavior as `simulateSolution`.

## `speval`

Evaluate state, input, dynamics, or input-derivative interpolation data from a processed solution.

Syntax:

```matlab
values = speval(solution, solType, nidx, T);
```

Inputs:

- `solution`: processed solution structure.
- `solType`: solution type. Common values include `'X'`, `'state'`, `'U'`, `'input'`, `'dX'`, `'dynamics'`, and `'dU'`.
- `nidx`: scalar or vector of state/input indices to evaluate.
- `T`: time vector where values should be evaluated.

Outputs:

- `values`: matrix with one column for each requested index and one row for each time value.

Side effects and assumptions:

- Requires interpolation fields such as `Xp`, `Up`, `dXp`, or `dUp` to exist in `solution`.
- Supports both standard piecewise-polynomial data and segmented LGR interpolation data.

## `genSolutionPlots`

Generate standard ICLOCS plots for a solved problem.

Syntax:

```matlab
genSolutionPlots(options, solution);
```

Inputs:

- `options`: options structure used for the solve.
- `solution`: solution returned by `solveMyProblem`.

Outputs:

- None.

Side effects and assumptions:

- Creates MATLAB figures for states, inputs, adjoints, error estimates, and constraint violations depending on the configured plot options.
- Supports both single-phase and multiphase solution structures.

## `checkProblem`

Validate the consistency of key problem-structure dimensions before solving.

Syntax:

```matlab
checkProblem(problem);
```

Inputs:

- `problem`: ICLOCS problem structure.

Outputs:

- None.

Side effects and assumptions:

- Throws an error when state, input, path-constraint, or boundary-constraint bounds and tolerances have inconsistent dimensions.
- Checks that the problem file follows the ICLOCS 2.5 data-field conventions.

## `setupIclocsPath`

Convenience helper for adding ICLOCS source folders to the MATLAB path.

Syntax:

```matlab
repoRoot = setupIclocsPath;
repoRoot = setupIclocsPath(repoRoot);
```

Inputs:

- `repoRoot`: optional path to the repository root.

Outputs:

- `repoRoot`: resolved repository root that was used to add paths.

Side effects and assumptions:

- Adds the repository `src` tree to the MATLAB path with `addpath(genpath(...))`.
- Does not add every example problem folder; add the specific example folder you want to run.
