function [ options, guess] = doWarmStart( options, guess, solution, data )
%doWarmStart - prepare the problem for warm start
%
% Syntax:  [ options, guess] = doWarmStart( options, guess, solution, data )
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------


if (isfield(options,'hpAdaptive') && options.hpAdaptive==1)
    options.start='Warm';
else
    options.start='Hot';
end


guess.states=solution.X;
guess.inputs=solution.U;
if isfield(solution,'Xp')
    guess.Xp=solution.Xp;
    guess.Up=solution.Up;
end


if isfield(solution,'TSeg_Bar')
    guess.TSeg_Bar=solution.TSeg_Bar;
end
guess.time=solution.T;
if strcmp(data.options.discretization,'hermite') || strcmp(data.options.discretization,'AutoDirect')
    guess.time_org=solution.T;
else
    guess.time_org=solution.coll.T;
end
guess.tf=solution.tf;
guess.t0=solution.t0;
if isfield(solution.multipliers,'lambdaNLP_1toN')
    guess.multipliers.lambda=solution.multipliers.lambdaNLP_1toN;
end
ng=deal(data.sizes{5});
nb=deal(data.sizes{6});
if ng
    if strcmp(data.options.discretization,'hermite') || strcmp(data.options.discretization,'AutoDirect')
        guess.timeFull=solution.coll.T;
        guess.multipliers.lambda_g=solution.multipliers.lambda_g;
    else
        guess.timeFull=solution.T;
        guess.multipliers.lambda_g=solution.multipliers.lambda_g;
    end
end
if nb
    guess.multipliers.lambda_b=solution.multipliers.lambda_b;
end
if isfield(solution.multipliers,'lambda_resconst')
    guess.multipliers.lambda_resconst=solution.multipliers.lambda_resconst;
end

if isfield(solution,'residuals') && strcmp(data.options.transcription,'integral_res_min')
    guess.residual=sqrt(solution.residuals.r);
end

if isfield(solution,'min_res_satisfy') && solution.min_res_satisfy==1
    guess.cost.ub=solution.cost;
end

end

