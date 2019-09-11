function bc=collmatch_endpt(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
%collmatch_endpt - for obtaining integrated residual minimization solution representation with matching of direct collocation result at endpoints
%
% Syntax:   bc=collmatch_endpt(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
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
X_endpoint=[x0';xf'];
x_endpoint=vdat.x_endpoint([1,end],:);
bc=X_endpoint(:)-x_endpoint(:);

end