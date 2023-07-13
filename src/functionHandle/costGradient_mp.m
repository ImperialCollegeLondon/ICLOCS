function grad_mp=costGradient_mp(z_mp,auxdata)
%COSTGRADIENT - Evaluate the gradient of the objective function, for a
%multi-phase problem
%
% Syntax:  grad_mp=costGradient_mp(z_mp,auxdata)
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

grad_mp=zeros(1,auxdata.mpdata.mpsizes.nz);

for i=1:length(auxdata.phasedata)
    data=auxdata.phasedata{i};
    z=zeros(data.nz,1);
    z(data.zidx.org.z,1)=z_mp(data.zidx.mp.z,1);
    
    switch data.options.transcription

        case {'multiple_shooting'}

            grad=multipleShooting('gradCost',z,data);

        % Direct Collocation
        case {'direct_collocation','direct_collocation_intres_reg'}

            switch data.options.discretization

                case{'discrete','euler','trapezoidal','hermite'}

                    grad=directCollocation('gradCost',z,data,i);

                case{'globalLGR','hpLGR'}

                    grad=directCollocationLGR('gradCost',z,data,i);

                case {'resMinConstructionForSolution'}

                    grad=resMinConstructionForSolution('gradCost',z,data);

                case {'resMinInterpolationForSolution'}

    %                 grad=resMinInterpolationForSolution('gradCost',z,data);
                    grad=resMinOptimization('gradCost',z,data);

            end

        case {'integral_res_min'}

            grad=resMinOptimization('gradCost',z,data);

    end
    
    grad_mp(1,data.zidx.mp.z)=grad_mp(1,data.zidx.mp.z)+grad(1,data.zidx.org.z);
    
end

end



%------------- END OF CODE --------------