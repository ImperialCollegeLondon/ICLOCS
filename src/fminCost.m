function [C GC HC]=fminCost(z,data)
%FMINCOST - Return the cost, the gradient and the hessian of the cost for fmincon
%
% Syntax:  [C GC HC] = fminCost(z,data)
%
% Inputs:
%    z     - Unknown NLP vector
%    data  - Data passed to the functions evaluated during optimization
%
% Outputs:
%    C     - Cost at point z
%    GC    - Gradient of the cost at z
%    HC    - Hessian of the cost at z
%
% Other m-files required: costFunction.m, costGradient.m, costHessian.m
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk


%------------- BEGIN CODE --------------


global sol

if nargout>1
    GC=costGradient(z,data);
    
    if nargout>2
        HC=costHessian(z,data,sol.JL);
    end
    
end


C=costFunction(z,data);