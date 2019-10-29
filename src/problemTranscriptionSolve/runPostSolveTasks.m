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
    solution.mp.cost=0;
    for i=1:length(problem.phases)
        [solution.phaseSol{i}]=postSolveAnalysis(problem.phases{i},solution.phaseSol{i},options.phaseoptions{i},data.phasedata{i});
    end
    solution.mp.cost=solution.mp.cost+solution.phaseSol{i}.cost;
    printSolveInfo(solution,options);
else
    if ~strcmp(options.discretization,'discrete')
        [solution]=postSolveAnalysis(problem,solution,options,data);
        printSolveInfo(solution,options);
    end
end

end

