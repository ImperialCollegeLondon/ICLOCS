 
function [pp,dpp]=Linsplines(T,X,F)
%  Linsplines - Approximate the trajectory in x with a linear
%             interpolation
%
% Syntax:  [pp,dpp]=Linsplines(T,X,F)
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

% Create the standard MATLAB piecewise polynomial (pp) structure for the
% linear spline
pp=struct('form','pp','breaks',T.','coefs',zeros(M,2),'pieces',M,'order',2,'dim',1);


% compute coefficients for the linear spline interpolation 
pp.coefs(:,2)=X(1:end-1); 
if length(X(1:end-1))==length(F(1:end-1))
	pp.coefs(:,1)=(X(2:end)-X(1:end-1))./DT;
    dpp=mkpp(T,F(1:end-1));
else
    pp.coefs(:,1)=(X(2:end)-X(1:end-1))./DT;
    dpp=[];
end

    
end
