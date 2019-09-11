function Y2 = legendreEval(A,X)
%legendreEval Evaluation Legendre polynomials
%
% Syntax:   Y2 = legendreEval(A,X)
%
% Inputs:
%    A - Matrix containing the fitted result from function legendrefit
%    X - The vector to be evalutated
%
% Outputs:
%    Y2 - The corresponding value after evaluation
%
% Adapted from LEGENDREFIT by
% Siqing Wu, <6sw21@queensu.ca> 
% Version: 1.2, Date: 2008-07-31
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


if size(X,2)~=1
    X=X';
end

N=length(A)-1;
%%% compute the Legendre polynomial coefficients matrix coeff
% coeff(i,j) gives the polynomial coefficient for term x^{j-1} in P_{i-1}(x)
if N > 1
    coeff = zeros(N+1);
    coeff([1 N+3]) = 1; % set coefficients of P_0(x) and P_1(x)
    % now compute for higher order: nP_n(x) = (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x)
    for ii = 3:N+1
        coeff(ii,:) = (2-1/(ii-1))*coeff(ii-1,[end 1:end-1]) - (1-1/(ii-1))*coeff(ii-2,:);
    end
else
    % simple case
    coeff = eye(N+1);
end

m = length(X);


%%% Evaluate the polynomials for every element in X
D = cumprod([ones(m,1) X(:,ones(1,N))], 2) * coeff.';

% % fitting (regression) result
Y2 = D*A; 

