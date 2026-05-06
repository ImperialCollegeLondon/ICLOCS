# Closed-Loop Simulation

This page summarizes the ICLOCS website guidance for closed-loop and NMPC simulations:

- [Closed-loop simulations](http://www.ee.ic.ac.uk/ICLOCS/GetStartedCLSim.html)

ICLOCS 2.5 supports closed-loop simulations with one or more MPC controller designs. The closed-loop examples live under `closedLoopExamples`.

## Controller Folders

Place each controller design in its own folder, for example:

- `MPCController1`
- `MPCController2`

Use unique file names inside each controller folder, commonly by appending `Phase1`, `Phase2`, and so on.

Each controller folder looks similar to a single-phase problem folder, except the main script is replaced by a controller script such as `ICLOCS_NLPSolver_Phase1.m`.

## Re-Optimization Triggers

The controller script decides when to solve again.

Standard update after the previous solution has become available:

```matlab
if simtime ~= 0 && simtime >= solution.t_ref + solution.computation_time
    re_opt = 1;
end
```

Event-based update:

```matlab
if simtime > solution.t_ref + 20
    re_opt = 1;
end
```

Constraint-based update can also be used when states approach configured bounds.

## Plant Dynamics

For Simulink simulations, plant dynamics can be supplied through an S-function script such as `Plant_sfun.m`. Configure:

- Number of continuous states.
- Number of outputs.
- Number of inputs.
- Plant equations in a function such as `Plant_eqn`.

## Initialization Script

The initialization script, typically `initSim.m`, adds controller folders to the MATLAB path, sets the simulation step, selects the initial mesh size, and fetches the first controller problem:

```matlab
addpath('MPCController1')
tstep = 0.1;
N_node = 20;
settings_func = @settings_myProblem_Phase1;
[problem,guess] = myProblem_Phase1;
```

## Controller Selection

The top-level `ICLOCS_NLPSolver.m` chooses which controller to call and configures:

- `opt_step`: interval between re-optimization attempts. Use `0` for update as soon as possible.
- `opt_t_min`: minimum horizon duration before falling back to open-loop implementation.

`opt_step` should be larger than the simulation time step when fixed-step re-optimization is used.

## Simulink Model

Update the Simulink model so the dynamics block points to the correct plant S-function. Run the simulation from the folder containing the initialization script and Simulink model.
