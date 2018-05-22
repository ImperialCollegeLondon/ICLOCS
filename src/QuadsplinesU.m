 
function [pp]=QuadsplinesU(T,U)

%  QuadsplinesU - Approximate the trajectory in u with a quadratic
%             interpolation
%
% Syntax:  [pp]=QuadsplinesU(T,U)
% 
% Inputs: 
%     T  - time vector   t=[t0, t1, ... ,tf]
%     U  - vector containing the input trajectory U=[u(t0),u(t1),...,u(tf)]
%
% Outputs:
%    pp  -   standard MATLAB piecewise polynomial (pp) structure containing the
%            quadratic interpolation of the discrete trajectory u  
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


DT=diff(T(1:2:end));      % differences  T(i+1)-T(i)
tk=zeros(size(T(1:2:end-1)));

M=(size(U,1)-1)/2;

uk=U(1:2:end-1);
ukh=U(2:2:end-1);
ukp1=U(3:2:end);

% Create the standard MATLAB piecewise polynomial (pp) structure for the
% quadratic spline
pp=struct('form','pp','breaks',T(1:2:end).','coefs',zeros(M,3),'pieces',M,'order',3,'dim',1);


% compute coefficients for the quadratic spline interpolation 
pp.coefs(:,3)=(2.*uk.*(DT + tk).*(DT./2 + tk))./DT.^2 - (4.*tk.*ukh.*(DT + tk))./DT.^2 + (2.*tk.*ukp1.*(DT./2 + tk))./DT.^2;
pp.coefs(:,2)=((4.*ukh.*(DT + tk))./DT.^2 - (2.*uk.*(DT + tk))./DT.^2 + (4.*tk.*ukh)./DT.^2 - (2.*tk.*ukp1)./DT.^2 - (2.*uk.*(DT./2 + tk))./DT.^2 - (2.*ukp1.*(DT./2 + tk))./DT.^2);
pp.coefs(:,1)=((2.*uk)./DT.^2 - (4.*ukh)./DT.^2 + (2.*ukp1)./DT.^2);



