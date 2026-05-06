# Options Reference

This is a compact repository copy of the options guidance from:

- [Single-phase options](http://www.ee.ic.ac.uk/ICLOCS/GetStartedOptionSingle.html)
- [Multi-phase options](http://www.ee.ic.ac.uk/ICLOCS/GetStartedOptionMultiPhase.html)

## Transcription

`options.transcription`:

- `direct_collocation`: direct collocation transcription.
- `integral_res_min`: integrated residual minimization, for supported single-phase workflows.

For integrated residual minimization, configure:

- `options.min_res_mode`: `alternating` or `weightedCost`.
- `options.min_res_priority`: `low_res_error` or `low_cost_value`.
- `options.resCostWeight`: penalty sequence for weighted-cost mode.
- `options.errortype`: `local_abs`, `int_res`, or `both`.

Multi-phase workflows currently use direct collocation.

## Discretization

`options.discretization`:

- `discrete`: discrete-time model.
- `euler`: Euler method.
- `trapezoidal`: trapezoidal method.
- `hermite`: Hermite-Simpson method.
- `globalLGR`: one global LGR polynomial.
- `hpLGR`: multiple LGR polynomial segments.
- `AutoDirect`: automatic direct-collocation choice.

## Result Representation

`options.resultRep`:

- `default`: representation matching the transcription method.
- `manual`: manually selected representation.
- `res_min`: representation by integrated residual minimization.
- `res_min_final_default`: residual-minimization representation only at the final mesh-refinement step.
- `res_min_final_manual`: final residual-minimization representation with manually selected intermediate representation.

Manual state representations include `linear`, `quadratic`, `cubic`, `pchip`, `Barycentric`, and `Legendre` depending on the discretization. Manual input representations include `constant`, `linear`, `quadratic`, `pchip`, `Barycentric`, and `Legendre` depending on the discretization.

## Derivatives

`options.derivatives`:

- `analytic`: use supplied analytic derivatives.
- `numeric`: finite-difference derivatives.
- `adigator`: algorithmic differentiation through ADiGator.

Analytic derivative mode requires the problem definition to link the relevant functions:

```matlab
problem.analyticDeriv.gradCost = @gradCost;
problem.analyticDeriv.hessianLagrangian = @hessianLagrangian;
problem.analyticDeriv.jacConst = @jacConst;
```

ADiGator mode requires `options.adigatorPath` to point to the directory containing `startupadigator.m`.

## NLP Solvers

`options.NLPsolver`:

- `ipopt`: recommended interior-point solver.
- `fmincon`: supported but usually slower; useful for sanity checks when Optimization Toolbox is available.
- `worhp`: supported for versions that still include a MATLAB interface.

Important IPOPT options include:

- `options.ipopt.tol`
- `options.ipopt.max_iter`
- `options.ipopt.hessian_approximation`

## Mesh Strategy

`options.meshstrategy`:

- `fixed`: fixed mesh.
- `mesh_refinement`: local mesh refinement iterations.
- `hp_flexible`: flexible adaptive mesh for hp-LGR direct collocation.

Mesh refinement choices include:

- `options.MRstrategy`: `aggressive` or `efficient`.
- `options.maxMRiter`.
- `options.mintimeinterval`.
- `options.maxtimeinterval`.
- `options.tau`: `0` for equispaced steps or a normalized nonuniform interval distribution.

## Other Settings

- `options.start`: `Cold`, `Warm`, or `Hot`.
- `options.scaling`: enable automatic scaling for variables with very different magnitudes.
- `options.resminEarlyStop`: allow early residual-minimization termination when tolerances are met.
- `options.ECH.enabled`: enable external constraint handling where supported.
- Regularization strategies include `off`, `reg_priority`, `MR_priority`, and `simultaneous` where configured.
