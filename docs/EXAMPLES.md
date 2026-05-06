# Examples

ICLOCS includes single-phase, multi-phase, fmincon-oriented, and closed-loop/NMPC examples. Many examples reuse the same function names within their own folders, so run one example folder at a time.

The full existing example matrix is available on the ICLOCS website at [www.ee.ic.ac.uk/ICLOCS/Examples.html](http://www.ee.ic.ac.uk/ICLOCS/Examples.html).

Expected qualitative outputs for selected examples are tracked in [`EXAMPLE_RESULTS.md`](EXAMPLE_RESULTS.md).

## Recommended Beginner Path

1. `exampleProblems/BangBang`
2. `exampleProblems/DoubleIntegratorTracking`
3. `exampleProblems/SecondOrderSingular`
4. `exampleProblems/CartPoleSwingUp`
5. `exampleProblems/GoddardRocket/SinglePhase`
6. `exampleProblems/GoddardRocket/MultiPhase`

## Example Matrix

| Example | Folder | Implementation | Solver | Notes |
| --- | --- | --- | --- | --- |
| Bang-bang Control | `exampleProblems/BangBang` | Easy | Easy | Detailed instructions; analytic derivatives; flexible hp-adaptive mesh demo. |
| Double Integrator Tracking | `exampleProblems/DoubleIntegratorTracking` | Easy | Easy | Detailed instructions; reference-tracking variant. |
| Hypersensitive Control | `exampleProblems/Hypersensitive` | Easy | Medium | Mesh-refinement demo. |
| Aly Chan Problem | `exampleProblems/Aly Chan` | Easy | Hard | Singular arc; residual minimization demo. |
| Van der Pol Oscillator | `exampleProblems/VanderPalOscillator` | Easy | Hard | Residual minimization demo. |
| Second Order Singular Control | `exampleProblems/SecondOrderSingular` | Easy | Medium | Rate-constraint demo. |
| Two-link Robot Arm | `exampleProblems/TwoLinkRobotArm` | Easy | Easy | Detailed instructions. |
| Cart Pole Swing-Up | `exampleProblems/CartPoleSwingUp` | Easy | Easy | Analytic derivatives; closed-loop simulation demo. |
| Batch Fermentor | `exampleProblems/BatchFermentor` | Easy | Hard | Detailed instructions; analytic derivatives. |
| Spaceship Control | `exampleProblems/SpaceshipControl` | Intermediate | Medium | Nonlinear dynamics. |
| Aircraft Minimum Fuel/Time Climb | `exampleProblems/SupersonicAircraftClimb` | Intermediate | Medium | Look-up tables and closed-loop simulation material. |
| Goddard Rocket, Single Phase | `exampleProblems/GoddardRocket/SinglePhase` | Intermediate | Hard | Residual minimization demo. |
| Orbit Raising | `exampleProblems/OrbitRaising` | Intermediate | Medium/Easy | Path and boundary constraints; alternative formulation demo. |
| Space Shuttle Reentry | `exampleProblems/SpaceShuttleReentry` | Hard | Medium | Detailed instructions. |
| Aircraft Go-around in Windshear | `exampleProblems/LandingWindshear` | Hard | Medium | Path and rate constraints. |
| Vehicle Parallel Parking | `exampleProblems/CarParking` | Hard | Hard | Automatic regularization demo. |
| Commercial Aircraft Flight Profile | `exampleProblems/CommericalFlightProfile` | Hard | Hard | External constraint handling demo. |
| Low-Thrust Orbit Transfer | `exampleProblems/LowThrustOrbitTransfer` | Hard | Hard | Long-horizon aerospace example. |
| Bang-bang Control, Multi-Phase | `exampleProblems/BangBang_TwoPhase` | Intermediate | Easy | Multi-phase formulation. |
| Goddard Rocket, Multi-Phase | `exampleProblems/GoddardRocket/MultiPhase` | Hard | Medium | Linkage and singular-arc condition material. |

## Running An Example

From the repository root in MATLAB:

```matlab
addpath('tools')
setupIclocsPath
exampleDir = fullfile(pwd,'exampleProblems','DoubleIntegratorTracking');
addpath(exampleDir)
cd(exampleDir)
main_DoubleIntegratorTracking
```

Each example usually contains:

- `main_*.m`: script that solves and plots the problem.
- `settings_*.m`: transcription, mesh refinement, derivative, and solver settings.
- `*_Dynamics_Internal.m`: dynamics used by the transcription.
- `*_Dynamics_Sim.m`: dynamics used for simulation and plotting.
- `gradCost_*.m`, `jacConst_*.m`, `hessianLagrangian_*.m`: optional analytic derivatives where provided.

## Adapting A Template

For a new single-phase problem, start from `usr/Single Phase Problem`. For a new multi-phase problem, start from `usr/Multi Phase Problem`.

Copy the template folder outside the repository or into a new working folder, then edit:

- `myProblem.m` or `myMultiPhaseProblem.m` for bounds, guesses, callbacks, and function handles.
- `settings_*.m` for solver, mesh, derivative, and output settings.
- `*_Dynamics_Internal*.m` for the optimal-control dynamics.
- `*_Dynamics_Sim*.m` for post-solve simulation.
- Analytic derivative files only if you enable analytic derivatives.

Leave the generated-function-handle section at the bottom of the problem file unchanged unless the template version changes.
