function data_org=checkDynamics( z0,data )
%checkDynamics - pre-testing of the user defined function, to catch errors more easily
%
% Syntax:  checkDynamics( z0,data )
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
% 
%------------- BEGIN CODE --------------

data_org=data;
if isfield(data,'dataNLP')
    data=data.dataNLP;
end
if (strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))
    data.sizes{5}=data.sizes{18}+data.sizes{19};
else  
    data.sizes{5}=data.sizes{15}+data.sizes{16};
end
data.options.transcription='direct_collocation';
data.options.derivatives='numeric';
data.data.resmin=0;
constraintFunction(z0,data);
nt=data.sizes{1};

z_zeros=zeros(size(z0));
data.QPSolverCheck= ~nt && (data.FD.FcnTypes.Ltype~=1) && (data.FD.FcnTypes.Etype~=1) && (data.FD.FcnTypes.Ftype~=1 && data.FD.FcnTypes.Ftype~=3) && (data.FD.FcnTypes.Gtype~=1 && data.FD.FcnTypes.Gtype~=3) && (data.FD.FcnTypes.Btype~=1 && data.FD.FcnTypes.Btype~=3);
if data.QPSolverCheck
    testout=typeTests(z_zeros,1,data);
    OSQP.P=sparse(size(testout.hessian.Lzz,1),size(testout.hessian.Lzz,2));
    OSQP.q=sparse(size(testout.gradCost.Lz,1),size(testout.gradCost.Lz,2));
    if data.FD.FcnTypes.Ltype==2 %linear cost
        OSQP.q=OSQP.q+testout.gradCost.Lz; %linear term 
    elseif data.FD.FcnTypes.Ltype==3 %quadratic cost
        OSQP.P=OSQP.P+testout.hessian.Lzz;
        OSQP.q=OSQP.q+testout.gradCost.Lz; %linear term
    end
    if data.FD.FcnTypes.Etype==2 %linear cost
        OSQP.P=OSQP.P+[]; %no quadratic term
        OSQP.q=OSQP.q+testout.gradCost.Ez; %linear term 
    elseif data.FD.FcnTypes.Etype==3 %quadratic cost
        OSQP.P=OSQP.P+testout.hessian.Lzz;
        OSQP.q=OSQP.q+testout.gradCost.Ez; %linear term
    end
    OSQP.A=[speye(length(z0));testout.jacConst.fzd;testout.jacConst.gzd;testout.jacConst.rcz;testout.jacConst.bz];
    ccst=[testout.const.fc;testout.const.gc;testout.const.avrcc;testout.const.bc];
    OSQP.l=[data.infoNLP.zl;data.infoNLP.cl-ccst]; %lower bound
    OSQP.u=[data.infoNLP.zu;data.infoNLP.cu-ccst]; %upper bound
end
data_org.OSQP=OSQP;


end

