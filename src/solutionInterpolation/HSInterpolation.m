 
function [pp,dpp]=HSInterpolation(T,x,dx)
%  Hsplines - Approximate the trajectory in x with a cubic Hermite
%             interpolation
%
% Syntax:  [pp,dpp]=Hsplines(T,x,dx)
% 
% Inputs: 
%     T  - time vector   t=[t0, t1, ... ,tf]
%     x  - vector containing the trajectory x=[x(t0),x(t1),...,x(tf)]
%     dx - vector containing the derivative of the trajectory dx/dt=[dx/dt(t0),dx/dt(t1),...,dx/dt(tf)]
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
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------


DT=diff(T(1:2:end));      % differences  T(i+1)-T(i)
% Dx=diff(x)./DT;  % differentiation (x(i+1)-x(i))/(T(i+1)-T(i))
M=size(DT,1);

X_k=x(1:2:end-1); 
F_k=dx(1:3:end-2);
F_kph=dx(2:3:end-1);
F_kp1=dx(3:3:end);

% Create the standard MATLAB piecewise polynomial (pp) structure for the
% cubic spline
pp=struct('form','pp','breaks',T(1:2:end).','coefs',zeros(M,4),'pieces',M,'order',4,'dim',1);


% compute coefficients for the cubic Hermite spline interpolation 
pp.coefs(:,4)=X_k;
pp.coefs(:,3)=F_k;
pp.coefs(:,2)=0.5*(-3*F_k+4*F_kph-F_kp1)./(DT);
pp.coefs(:,1)=1/3*(2*F_k-4*F_kph+2*F_kp1)./(DT.*DT);


% Create the standard MATLAB piecewise polynomial (pp) containing the estimation
% of the derivatives of x

dpp=pp;
dpp.coefs=repmat([3,2,1],pp.pieces,1).*pp.coefs(:,1:end-1);
dpp.order=pp.order-1;