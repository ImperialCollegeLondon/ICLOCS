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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------


switch data.options.transcription
    
    case {'multiple_shooting'}
        
        grad=multipleShooting('gradCost',z,data);
        
    % Direct Collocation
    case{'discrete','euler','trapezoidal','hermite'}
        
        grad=directCollocation('gradCost',z,data);
        
    case{'globalLGR','hpLGR'}
        
        grad=directCollocationLGR('gradCost',z,data);

end

%------------- END OF CODE --------------