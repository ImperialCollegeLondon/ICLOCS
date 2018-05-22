 
function [pp,dpp]=Quadsplines(T,X,F)

%  Quadsplines - Approximate the trajectory in x with a quadratic
%             interpolation
%
% Syntax:  [pp,dpp]=QuadsplinesU(T,X,F)
% 
% Inputs: 
%     T  - time vector   t=[t0, t1, ... ,tf]
%     X  - vector containing the trajectory X=[X(t0),X(t1),...,X(tf)]
%     F - vector containing the derivative of the trajectory F=[f(t0),f(t1),...,f(tf)]
%
% Outputs:
%    pp  -   standard MATLAB piecewise polynomial (pp) structure containing the
%            Hermite interpolation of the discrete trajectory x   
%    dpp -   standard MATLAB piecewise polynomial (pp) structure containing the
%            the derivative of the Hermite interpolation of x  
%
% Other m-files required: none
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


  
DT=diff(T);      % differences  T(i+1)-T(i)

M=size(X,1)-1;

xk=X(1:end-1);
fk=F(1:end-1);
fkp1=F(2:1:end);

% Create the standard MATLAB piecewise polynomial (pp) structure for the
% quadratic spline
pp=struct('form','pp','breaks',T.','coefs',zeros(M,3),'pieces',M,'order',3,'dim',1);


% compute coefficients for the quadratic spline interpolation 
pp.coefs(:,3)=xk; 
pp.coefs(:,2)=fk;
pp.coefs(:,1)=0.5*(-fk+fkp1)./DT;

% Create the standard MATLAB piecewise polynomial (pp) containing the estimation
% of the derivatives of x

dpp=pp;
dpp.coefs=repmat([2,1],pp.pieces,1).*pp.coefs(:,1:end-1);
dpp.order=pp.order-1;
