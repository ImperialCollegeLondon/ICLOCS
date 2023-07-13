function jacConst=constraintJacobian(z,data)
%CONSTRAINTJACOBIAN - Evaluate the Jacobian of the constraints
%
% Syntax: jacConst=constraintJacobian(z,data)
%
% Inputs:
%    z - NLP variable 
%    data - Structure with data needed to evaluate CostFn
%
% Outputs:
%    jacConst - Sparse matrix with jacobian of the constraints at z
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
        
         jacConst=multipleShooting('jacConst',z,data);
   
    case {'direct_collocation','direct_collocation_intres_reg'}
        
        switch data.options.discretization
            
           case{'globalLGR','hpLGR'}
               
                 jacConst=directCollocationLGR('jacConst',z,data,1);

            case {'resMinConstructionForSolution'}

                jacConst=resMinConstructionForSolution('jacConst',z,data);

            case {'resMinInterpolationForSolution'}

%                 jacConst=resMinInterpolationForSolution('jacConst',z,data);
                jacConst=resMinOptimization('jacConst',z,data,1);


            otherwise % Direct Transcription Method

                 jacConst=directCollocation('jacConst',z,data,1);
         
        end
        
    case {'integral_res_min'}
        
        jacConst=resMinOptimization('jacConst',z,data,1);

end

 %------------- END OF CODE --------------                     






