function b=callback_minres(t,f,vdat)
%callback - customized callback function for integral residual minimization
%
% Syntax:  b=callback_minres(t,f,vdat)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

global sol

solution=sol{1};

% check prime feasible of the problem
flag_feas=0;
if all(solution.const<=(vdat.infoNLP.cu+vdat.infoNLP.c_Tol)) && all(solution.const>=(vdat.infoNLP.cl-vdat.infoNLP.c_Tol)) && all(solution.z<=vdat.infoNLP.zu) && all(solution.z>=vdat.infoNLP.zl)
    solution.lastfeasible=solution.z;
    flag_feas=1;
end
        
% Check residual error tolerance
if all(solution.constRes<vdat.dataNLP.data.discErrorTol_FullScaling)
    if flag_feas
        b = false;
    else
        b = true;
    end
else
    b = true;
end
