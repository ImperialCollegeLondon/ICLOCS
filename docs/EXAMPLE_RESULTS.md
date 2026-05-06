# Example Results

This page captures expected qualitative outputs from the ICLOCS website so users know what a successful run should look like before they compare exact numerical values. Computation times are historical reference values from the website and depend strongly on machine, MATLAB version, solver build, derivative mode, and tolerances.

## Bang-Bang Control

Folder: `exampleProblems/BangBang`

Expected behavior:

- State trajectory moves from the initial position to the terminal position.
- Control is bang-bang and switches near the middle of the trajectory.
- Mesh refinement places additional mesh points near the switching structure.

Website reference:

- Hermite-Simpson after one mesh-refinement iteration.
- About 0.7 seconds with analytic derivatives on the reference desktop.
- Flexible hp-adaptive mesh captures the control switch near 20 seconds.

Source: [Bang-bang Control](http://www.ee.ic.ac.uk/ICLOCS/ExampleBangBang.html)

## Orbit Raising

Folder: `exampleProblems/OrbitRaising`

Expected behavior:

- Radius increases toward the final orbit.
- Open-loop simulation tracks the optimized state and input trajectories.
- Alternative angle-control formulation should be faster than the original direction-vector formulation.

Website reference:

- Hermite-Simpson with mesh refinement starting from 40 mesh points.
- About 12 seconds for the original formulation with finite-difference derivatives on the reference desktop.
- About 1.4 seconds for the alternative formulation with analytic derivatives.
- Nondimensional time corresponds to approximately 193 days.

Source: [Orbit Raising](http://www.ee.ic.ac.uk/ICLOCS/ExampleOrbitRaising.html)

## Low-Thrust Orbit Transfer

Folder: `exampleProblems/LowThrustOrbitTransfer`

Expected behavior:

- Long-horizon trajectory with several state and input plots.
- Open-loop simulation is shown against the optimized solution.

Website reference:

- Hermite-Simpson with mesh refinement starting from 150 mesh points.
- About 373 seconds after one mesh-refinement iteration with finite-difference derivatives on the reference desktop.

Source: [Low-Thrust Orbit Transfer](http://www.ee.ic.ac.uk/ICLOCS/ExampleLowThrustOrbitTransfer.html)

## Using These References

Use these notes as smoke-test expectations, not strict regression baselines. Future regression tests should store numerical tolerances separately and should be gated on solver availability.
