function [dL,dE]=gradCost_BatchFermentor(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,t0,tf,data)

%GRADCOST - Return the gradient of the cost in analytic form
%
% Syntax:  [dL,dE]=gradCost(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    see directCollocation.m
%
% Outputs:
%   dL - Gradient of the stage cost L  wrt. t,p,x,u 
%   dE - Gradient of E with respect to tf,p,x0,u0,uf,xf
%  
%  Gradient dL:
%  If the analytic form for L is supplied, set  dL.flag=1 otherwise set
%  dL.flag=0. If tf is a variable of the problem  and dL.flag=1 it is necessary to
%  specify  the derivative of L with respect to time t. 
%  You need to specify dL.dx and dL.du whenever dL.flag=1. 
%  dL.dx, dL.du,  dL.dp and dL.dt must be vectorized.
%  The derivative of L with respect to the i-th state variable, evaluated along all the horizon,
%  corresponds to the i-th column of dL.dx. The same rule holds for dL.du, dL.dp, dL.dt 
%  The i-th state and input
%  xi and ui, evaluated at the time instants t=[t0,...tk,...tf], are column vectors 
%  taken as X(:,i) and U(:,i)
%  The i-th state and input references  xri and uri for the i-th variable are column vectors taken 
%   as Xr(:,i) and Ur(:,i)  
%    
%  Gradient dE:
%  If the analytic form for E is supplied set  dE.flag=1 otherwise set
%  dE.flag=0. If the derivative of E does not depend on some variable it is possible to set
%  the respective variable as the empty matrix. For instance if it does not depend on x0 set
%  dE.dx0=[];
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------
[ X,Xr,U,Ur,P,x0,xf,u0,uf,p,data ] = batchScaleBack(X,Xr,U,Ur,P,x0,xf,u0,uf,p,data);

% Stage cost gradients
n=deal(data.sizes{3});
Lt=ones(size(t));

dL.flag=1; 
dL.dp=[];
dL.dt=[0*t];
dL.dx=kron(zeros(1,n),Lt);
dL.du=2*0.00001*U;



% Terminal cost gradients
dE.flag=1;
dE.dt0=[];
dE.dtf=[];
dE.dp=[];
dE.dx0=[];
dE.du0=[];
dE.dxf=[0 -xf(4) 0 -xf(2)];
dE.duf=[];

%------------- END CODE --------------



