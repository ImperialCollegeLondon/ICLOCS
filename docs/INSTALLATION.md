# Installation

These instructions are for users installing ICLOCS from a fresh clone or download.

The existing ICLOCS website installation page is available at [www.ee.ic.ac.uk/ICLOCS/Downloads.html](http://www.ee.ic.ac.uk/ICLOCS/Downloads.html).

## Requirements

- MATLAB R2017a or later.
- A nonlinear programming solver supported by ICLOCS. Most examples default to IPOPT through `ipopt.mex`.
- Optional: Optimization Toolbox if you want to use `fmincon`.
- Optional: ADiGator if you want algorithmic differentiation workflows.
- Optional: WORHP if you have a supported version with the MATLAB interface.

Because this repository has existed for a long time and examples vary, check the `settings_*.m` file for the example you plan to run. That file shows the selected derivative method, solver, mesh strategy, and solver-specific options.

## Clone Or Download

Clone the repository:

```bash
git clone https://github.com/ImperialCollegeLondon/ICLOCS.git
cd ICLOCS
```

Or download the repository ZIP from GitHub and extract it.

## Add ICLOCS To The MATLAB Path

From the repository root in MATLAB:

```matlab
addpath('tools')
setupIclocsPath
```

`setupIclocsPath` adds the ICLOCS `src` tree to the MATLAB path. Add only the example folder you are working on to avoid path conflicts between examples that intentionally reuse function names:

```matlab
addpath(fullfile(pwd,'exampleProblems','DoubleIntegratorTracking'))
```

## Configure A Solver

IPOPT is the preferred solver for most ICLOCS workflows. It is open source, but MATLAB users still need a working MATLAB interface such as `ipopt.mex`. Windows users have historically used OPTI Toolbox as one route to obtaining IPOPT.

`fmincon` is supported but not recommended for performance-sensitive solves. It can be useful as a sanity check when Optimization Toolbox is installed:

```matlab
options.NLPsolver = 'fmincon';
```

WORHP is also supported for versions that include the MATLAB interface. Newer WORHP releases may not include that interface.

## Run A First Example

Start with the Double Integrator Tracking example:

```matlab
addpath('tools')
setupIclocsPath
exampleDir = fullfile(pwd,'exampleProblems','DoubleIntegratorTracking');
addpath(exampleDir)
cd(exampleDir)
main_DoubleIntegratorTracking
```

If IPOPT is not installed, open `settings_DoubleIntegratorTracking.m` and review the `options.NLPsolver` section before running the example.

## Troubleshooting

- `Undefined function 'solveMyProblem'`: run `setupIclocsPath` from the repository root.
- `Undefined function 'ipopt'`: install IPOPT for MATLAB or configure the example to use a solver available on your machine.
- `Undefined function 'optimoptions'`: the selected settings file references Optimization Toolbox.
- Unexpected functions are called: remove broad example paths and add only the one example directory you are running.
- ADiGator errors: check `options.derivatives` and `options.adigatorPath` in the example settings file.
