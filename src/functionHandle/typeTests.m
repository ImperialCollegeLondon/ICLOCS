function testout=typeTests(z,sigma,data)
%COMPUTEHESSIAN - Evaluate the Hessian of the Lagrangian at z
%
% Syntax: hessian=computeHessian(z,sigma,lambda,data)
%
% Inputs:
%    z - NLP variable 
%    sigma - Cost function factor
%    lambda - Vector of Lagrange multipliers
%    data - Structure with data needed to evaluate CostFn
%
% Outputs:
%    hessian - Sparse lower triangular segment of the hessian at z 
%
% Other m-files required: multipleShooting.m or directCollocation.m
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

data.sigma=sigma;

switch data.options.transcription
    
    case {'direct_collocation','direct_collocation_intres_reg'}
        
        switch data.options.discretization
            
            case{'globalLGR','hpLGR'} % p/hp Transcription Method

                testout.gradCost=directCollocationLGR('gradCostTest',z,data,1);
                testout.hessian=directCollocationLGR('hessianTest',z,data,1);
                testout.jacConst=directCollocationLGR('jacConstTest',z,data,1);
                testout.const=directCollocationLGR('constTest',z,data,1);
            otherwise % h Transcription Method
                testout.gradCost=directCollocation('gradCostTest',z,data,1);
                testout.hessian=directCollocation('hessianTest',z,data,1);
                testout.jacConst=directCollocation('jacConstTest',z,data,1);
                testout.const=directCollocation('constTest',z,data,1);
        end
          

end

 %------------- END OF CODE --------------                     
