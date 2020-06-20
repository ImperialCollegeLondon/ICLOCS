function [solution] = runPostSolveTasks(problem,solution,options,data)
%runPostSolveTasks - run Post-Solve Tasks
%
% Syntax:  [solution] = runPostSolveTasks(problem,solution,options,data)
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


if isfield(solution,'mp')
    solution.mp.cost.J=0;
    for i=1:length(problem.phases)
        if solution.phaseSol{i}.tf-solution.phaseSol{i}.t0==0
            warning('Certain phases became redundant (i.e. tf-t0=0), post processing of solution might fail');
        else
            [solution.phaseSol{i}]=postSolveAnalysis(problem.phases{i},solution.phaseSol{i},options.phaseoptions{i},data.phasedata{i});
            solution.mp.cost.J=solution.mp.cost.J+solution.phaseSol{i}.cost.J;
        end
    end
    printSolveInfo(solution,options);
else
    [solution]=postSolveAnalysis(problem,solution,options,data);
    if ~strcmp(options.discretization,'discrete')
        printSolveInfo(solution,options);
    end
end

end

