 
function [Up,dUp]=HSInterpolationU(T,u)
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

U_k=u(1:2:end-1);
U_kph=u(2:2:end-1);
U_kp1=u(3:2:end);

% Create the standard MATLAB piecewise polynomial (pp) containing the estimation
% of the derivatives of x

Up=struct('form','pp','breaks',T(1:2:end).','coefs',zeros(M,3),'pieces',M,'order',3,'dim',1);

% compute coefficients for the cubic Hermite spline interpolation 
Up.coefs(:,3)=U_k;
Up.coefs(:,2)=(-3*U_k+4*U_kph-U_kp1)./DT;
Up.coefs(:,1)=(2*U_k-4*U_kph+2*U_kp1)./(DT.*DT);

dUp=struct('form','pp','breaks',T(1:2:end).','coefs',zeros(M,2),'pieces',M,'order',2,'dim',1);
% compute coefficients for the cubic Hermite spline interpolation 
dUp.coefs(:,2)=(-3*U_k+4*U_kph-U_kp1)./(DT);
dUp.coefs(:,1)=(2*U_k-4*U_kph+2*U_kp1)./(DT);
