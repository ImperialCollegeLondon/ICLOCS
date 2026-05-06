# Problem Formulation

This guide condenses the existing ICLOCS website guidance into repository documentation. The original pages are still useful references:

- [Single-phase problem formulation](http://www.ee.ic.ac.uk/ICLOCS/GetStartedProblemSingle.html)
- [Multi-phase problem formulation](http://www.ee.ic.ac.uk/ICLOCS/GetStartedProblemMultiPhase.html)

## Single-Phase Structure

Start from `usr/Single Phase Problem` when creating a new single-phase optimal-control problem. A typical problem folder contains:

- `myProblem.m`: problem bounds, guesses, cost functions, boundary constraints, and function handles.
- `myProblem_Dynamics_Internal.m`: dynamics and path constraints used by transcription.
- `myProblem_Dynamics_Sim.m`: optional simulation dynamics.
- `settings_myProblem.m`: transcription, mesh, derivative, and solver settings.
- Optional analytic derivative files: `gradCost_myProblem.m`, `jacConst_myProblem.m`, and `hessianLagrangian_myProblem.m`.

## Internal Dynamics

The internal dynamics function receives state, input, static parameter, time, and data arguments:

```matlab
function [dx,g_eq,g_neq] = myProblem_Dynamics_Internal(x,u,p,t,vdat)
```

Use vectorized expressions where possible. Return:

- `dx`: right-hand side of the ODE.
- `g_eq`: equality path constraints, optional.
- `g_neq`: inequality path constraints, optional.

## Problem Definition

In `myProblem.m`, define:

- Function handles for internal dynamics and simulation dynamics.
- Optional analytic derivative function handles.
- Initial and final time bounds.
- Static parameter bounds and guesses.
- Initial, path, and terminal state bounds.
- State error and box-constraint tolerances.
- Control/input bounds, first-control bounds, optional rate bounds, and tolerances.
- Path constraint counts, bounds, and tolerances.
- Boundary constraint bounds and tolerances.
- Auxiliary data in `problem.data`.
- Lagrange/stage cost in `L_unscaled`.
- Mayer/boundary cost in `E_unscaled`.
- Additional boundary constraints in `b_unscaled`.

The template section that assembles function handles and returns to `main.m` should be left unchanged unless the template itself is updated.

## Multi-Phase Structure

Start from `usr/Multi Phase Problem` for linked phases. A multi-phase setup has:

- `myMultiPhaseProblem.m`: shared time variables, shared parameters, linkage constraints, and phase initialization.
- One problem definition file per phase.
- One settings file per phase.
- A multi-phase settings file for solver-wide choices.

Each phase follows the same broad pattern as a single-phase problem, except phase start and end times refer to shared indices in `problem.mp.time`.

## Linkage Constraints

Multi-phase problems connect phases through linkage constraints. ICLOCS 2.5 separates linear and nonlinear linkage constraints:

```matlab
problem.mp.constraints.bll.linear = [...];
problem.mp.constraints.blu.linear = [...];
problem.mp.constraints.blTol.linear = [...];

problem.mp.constraints.bll.nonlinear = [...];
problem.mp.constraints.blu.nonlinear = [...];
problem.mp.constraints.blTol.nonlinear = [...];
```

The equations are defined in `bclink(x0,xf,u0,uf,p,t0,tf,vdat)`. Phase values are accessed through cells, for example `xf{4}(2)` for the second terminal state in phase 4.

## Solving

Single-phase:

```matlab
[problem,guess] = myProblem;
options = problem.settings(30);
[solution,MRHistory] = solveMyProblem(problem,guess,options);
```

Multi-phase:

```matlab
[problem,guess,options.phaseoptions] = myMultiPhaseProblem;
options.mp = settings_myMultiPhaseProblem;
[solution,MRHistory] = solveMyProblem(problem,guess,options);
```
