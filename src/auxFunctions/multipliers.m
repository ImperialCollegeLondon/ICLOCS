function [lambdai,Time]=multipliers(solution,data)
 
% MULTIPLIERS - Extract the adjoint variables and the multipliers relative to
%               the path constraints and the boundary constraints from the
%               optimal solution
%
%
% Syntax:  [lambdai,Time]=multipliers(solutions,data)
%
% Inputs:
%    solutions - Structure containing the solution of the optimization problem
%    data - structure containing matrices to format data 
%
% Outputs: 
%    lambdai - structured variable containing multipliers
%
%              lambdai.x0: multiplier of the constraint on (data.x0t-data.x0) 
%                          it is introduced in the Lagrangian when
%                          data.cx0=1 (see user guide)
%              lambdai.adjoint: adjoint variables 
%              lambdai.g: multipliers associated with the path constraints
%                         and they are not present if there are not path
%                         constraints (ng=0)
%              lambdai.b: multipliers associated with the boundary
%                         constraints if presents
%     Time -  structured variable containing the time instants relative  to
%             the multipliers
%
%             Time.adjoint: time instants  where the adjoint variables are evaluated. 
%             Time.g:       time instants  where the path constraints are evaluated if presents. 
%
% Other m-files required:  standard
% MAT-files required: none
%
% Copyright (C) 2013 Paola Falugi, Eric Kerrigan. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 15 May 2013
% iclocs@imperial.ac.uk



% Define some variables
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{:});
vdat=data.data;
t0=data.t0;
ngb=n*M+ng*M;
z=solution.z;


if nt; tf=z(1);else   tf=data.tf(1);end

T=(tf-t0)*[0;cumsum(data.tau)]*data.Nm/ns+data.k0;
Time.adjoint=T(ns)+T(1:ns:end-1);



switch(data.options.NLPsolver)
case{'ipopt'} 
  lambda_midpoint=reshape(solution.multipliers.lambda(n+1:n*M),n,M-1)';
case{'fmincon'}  
  lambda_midpoint=reshape(solution.multipliers.eqnonlin(n+1:n*M),n,M-1)';
end    
  lambda_midpoint=lambda_midpoint(ns:ns:M-1,:);
% if M<=2
%       lambda=lambda_midpoint;
%  else
%       lambda=interp1((1.5:(M-1)/ns+1)',lambda_midpoint,(1:(M-1)/ns+1)','linear','extrap');
% end 

lambdai.x0=solution.multipliers.lambda(1:n)';
lambdai.adjoint=lambda_midpoint;
switch(data.options.NLPsolver)
case{'ipopt'} 
if ng
  lambdai.g=reshape(solution.multipliers.lambda(n*M+1:ngb),ng,M)';  
  Time.g=T(1:ns:end);
end

if nb
  lambdai.b=solution.multipliers.lambda(ngb+1:ngb+nb); 
end  
end    