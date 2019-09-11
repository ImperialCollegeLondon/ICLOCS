function grad=costGradient(z,data)
%COSTGRADIENT - Evaluate the gradient of the objective function
%
% Syntax:  grad = costGradient(z,data)
%
% Inputs:
%    z    - NLP variable 
%    data - Structure with data needed to evaluate GradCostFn
%
% Outputs:
%    grad - Vector of values of the gradient of the cost at z
%
% Other m-files required: multipleShooting.m or directTranscription.m
% Subfunctions: none
% MAT-files required: none
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


switch data.options.transcription
    
    case {'multiple_shooting'}
        
        grad=multipleShooting('gradCost',z,data);
        
    % Direct Collocation
    case {'direct_collocation'}
        
        switch data.options.discretization

            case{'discrete','euler','trapezoidal','hermite'}

                grad=directCollocation('gradCost',z,data,1);

            case{'globalLGR','hpLGR'}

                grad=directCollocationLGR('gradCost',z,data,1);

            case {'resMinConstructionForSolution'}

                grad=resMinConstructionForSolution('gradCost',z,data);

            case {'resMinInterpolationForSolution'}

%                 grad=resMinInterpolationForSolution('gradCost',z,data);
                grad=resMinOptimization('gradCost',z,data,1);
        
        end
        
    case {'integral_res_min'}
        
        grad=resMinOptimization('gradCost',z,data,1);

end

%------------- END OF CODE --------------