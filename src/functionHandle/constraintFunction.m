function constraints=constraintFunction(z,data)
%CONSTRAINTFUNCTION - Evaluate the constraint functions
%
% Syntax:  constraints = constraintFunction(z,data)
%
% Inputs:
%    z    - NLP variable 
%    data - Structure with data needed to evaluate ConstFn
%
% Outputs:
%    constraints - Vector of values of the constraint function at x
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
        
        constraints=multipleShooting('const',z,data);
    
    % Direct Collocation
    case {'direct_collocation','direct_collocation_intres_reg'}
        
        switch data.options.discretization
            
            case{'discrete','euler','trapezoidal','hermite'}

                constraints=directCollocation('const',z,data,1);

            case{'globalLGR','hpLGR'}

                constraints=directCollocationLGR('const',z,data,1);


            case {'resMinConstructionForSolution'}

                constraints=resMinConstructionForSolution('const',z,data);

            case {'resMinInterpolationForSolution'}

                constraints=resMinOptimization('const',z,data,1);
        end
        
    case {'integral_res_min'}
        
        constraints=resMinOptimization('const',z,data,1);
        
end

%------------- END OF CODE --------------