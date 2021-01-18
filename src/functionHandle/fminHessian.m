function [Hout]=fminHessian(z,lambda,data,nlp)
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
n_neq_lb=length(nlp.ind_ineq_lb);
n_neq_ub=length(nlp.ind_ineq_ub);
lambda_full=zeros(nlp.nnod+length(nlp.ind_eq)+n_neq_lb+n_neq_ub,1);
lambda_full([nlp.ind_eqODE,nlp.nnod+nlp.ind_eq'])=lambda.eqnonlin;
lambda_full(nlp.nnod+nlp.ind_ineq_ub)=lambda_full(nlp.nnod+nlp.ind_ineq_ub)+lambda.ineqnonlin(1:n_neq_ub);
lambda_full(nlp.nnod+nlp.ind_ineq_lb)=lambda_full(nlp.nnod+nlp.ind_ineq_lb)-lambda.ineqnonlin(1+n_neq_ub:n_neq_ub+n_neq_lb);
 
Hes=computeHessian(z,1,lambda_full,data);          % Evaluate the constraints
Hout = triu(Hes.',1) + tril(Hes) ;



