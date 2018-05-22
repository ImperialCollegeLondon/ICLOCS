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
        
        constraints=multipleShooting('const',z,data);
    
    % Direct Collocation
    case{'discrete','euler','trapezoidal','hermite'}
        
        constraints=directCollocation('const',z,data);
        
    case{'globalLGR','hpLGR'}
        
        constraints=directCollocationLGR('const',z,data);

end

%------------- END OF CODE --------------