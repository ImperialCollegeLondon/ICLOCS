function [  Xout] = speval( Xp,n,varargin)
%speval - evaluate the spline functions
%
% Syntax:   Xout = speval( Xp,n,varargin)
%
% Inputs:
%    Xp  - Structure containing polynomial coefficients.
%    n - Return result for the nth state/input variable
%    varargin congtains
%    TSeg_Bar - time interval barrier (p/hp only)
%    T - Time vector
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

if nargin==4
    TSeg_Bar=varargin{1};
    T=varargin{2};
else
    T=varargin{1};
end

if isstruct(Xp{1})
    Xout=ppval(Xp{n},T);
else
    Xout=zeros(length(T),1);Xpi=Xp(:,n);
    for i=1:length(TSeg_Bar)-1 
        Tau_Seg=normalizeT(T(T>=TSeg_Bar(i) & T<TSeg_Bar(i+1)),TSeg_Bar(i),TSeg_Bar(i+1));
        Xout((T>=TSeg_Bar(i) & T<TSeg_Bar(i+1)),:)=legendreEval(Xpi{i},Tau_Seg);
    end
    if TSeg_Bar(end)==T(end)
        Xout(end)=legendreEval(Xpi{i},1);
    end
end

end

