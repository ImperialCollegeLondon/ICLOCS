function Y2 = BaryEval(A,x)
%BaryEval Barycentric Lagrange Interpolation.
%
% Syntax:   Y2 = BaryEval(A,X)
%
% Inputs:
%    A - First column the LGR points, Second column the data
%    X - The vector to be evalutated
%
% Outputs:
%    Y2 - The corresponding value after evaluation
%
% Adapted from barylag.m by Greg von Winckel
% Original Copyright info:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% barylag.m
%
% Interpolates the given data using the Barycentric
% Lagrange Interpolation formula. Vectorized to remove all loops
%
% data - a two column vector where column one contains the
%        nodes and column two contains the function value 
%        at the nodes
% p - interpolated data. Column one is just the 
%     fine mesh x, and column two is interpolated data 
%
% Reference:
%
% (1) Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange 
%     Interpolation" 
%     http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/berrut.ps.gz
% (2) Walter Gaustschi, "Numerical Analysis, An Introduction" (1997) pp. 94-95
%
%
% Written by: Greg von Winckel       03/07/04
% Contact:    gregvw@chtm.unm.edu 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

if ~isempty(x)
    M=size(A,1);     N=length(x);

    if isrow(x)
        x=x';
    end

    % Compute the barycentric weights
    X=repmat(A(:,1),1,M);

    % matrix of weights
    W=repmat(1./prod(X-X.'+eye(M),1),N,1);

    % Get distances between nodes and interpolation points
    xdist=repmat(x,1,M)-repmat(A(:,1).',N,1);

    % Find all of the elements where the interpolation point is on a node
    [fixi,fixj]=find(xdist==0);

    % Use NaNs as a place-holder
    xdist(fixi,fixj)=NaN;
    H=W./xdist;

    % Compute the interpolated polynomial
    Y2=(H*A(:,2))./sum(H,2);

    % Replace NaNs with the given exact values. 
    Y2(fixi)=A(fixj,2);
else
    Y2=[];
end

