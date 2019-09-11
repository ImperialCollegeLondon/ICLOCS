function [C,Ceq,jacC,jacCeq]=fminConst(z,data,nlp)
%FMINCONST - Return the constraints and jacobian of the constraints for fmincon
%
% Syntax:  [C Ceq jacC jacCeq]=fminConst(z,data,NLP)
%
% Inputs:
%    z     - Unknown NLP vector
%    data  - Data passed to the functions evaluated during optimization
%    nlp   - Data passed to the functions evaluated during optimization
%
% Outputs:
%    C       - Inequality constraints
%    Ceq     - Equality constraints
%    jacC    - Jacobian of the inequality constraints
%    jacCeq  - Jacobian of the equality constraints
%
% Other m-files required: costFunction.m, costGradient.m
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



[nt,np,n,m,ng,nb,M]=deal(data.sizes{1:7}); % Get some parameters
const=constraintFunction(z,data);          % Evaluate the constraints
jac=constraintJacobian(z,data);          % Evaluate the constraint jacobian
nnod=n*M;
nci=(ng || nb);

Ceq=const(1:nnod);

if nci
 flag_eq=(nlp.cu(nnod+1:end)-nlp.cl(nnod+1:end)==0);
 ind_eq=find(flag_eq); 
 ind_ineq=find(~flag_eq);
 Ceq=[Ceq;const(nnod+ind_eq)-nlp.cu(nnod+ind_eq)]; 
 C=[const(nnod+ind_ineq)-nlp.cu(nnod+ind_ineq);-const(nnod+ind_ineq)+nlp.cl(nnod+ind_ineq)];
else
 C=[];   
end    
 
 


if nargout>2
    jacCeq=jac(1:nnod,:)';
    if  nci
     jacCeq=[jacCeq,jac(nnod+ind_eq,:)' ];
     jacC=jac(nnod+ind_ineq,:)'; 
    else
     jacC=[];
    end
end


