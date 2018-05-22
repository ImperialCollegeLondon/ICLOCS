function varargout = legendrefit(X, Y, N, method)
%legendrefit Legendre polynomial Fitting
%
% Syntax:   Y2 = legendreEval(A,X)
%
% Inputs:
%    X, Y - Data for fitting
%    N - Polynomial order used for fitting
%    method - method for finding the solution in least squares sense
%       * inv - via direct inversion (default)
%       * qr - via QR decomposition
%       * chol - via Cholesky decomposition
%
% Outputs:
%    varargout - The fitted results
%
% Adapted from LEGENDREFIT by
% Siqing Wu, <6sw21@queensu.ca> 
% Version: 1.2, Date: 2008-07-31
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk


% make sure Y is a vector
if all(size(Y)>1)
    error('Input data should be a vector.')
end

% check N
if nargin < 2
    fprintf('Order N is not specified. Default order N=2 is used.\n')
    N = 2; % set N to default
end
if (N < 0 || round(N) ~= N)
    fprintf('Input N=%g is not valid. Default order N=2 is used.\n', N)
    N = 2;
end

% default method for finding the solution in least squares sense
if nargin < 3
    method = 'inv';
end

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

m = length(Y);
Y = Y(:); % make it a column


%%% Evaluate the polynomials for every element in X
D = cumprod([ones(m,1) X(:,ones(1,N))], 2) * coeff.';


%%% Find weighting coefficients for the linear combination of polynomials
switch method
    case 'inv'
        % Solution 1 (default)
        % for ill-conditioned case, method 'qr' may be a better choice
        A = (D.'*D)\(D.'*Y); % inverting the normal equations matrix directly 

    case 'chol'
        % Solution 2: Cholesky decomposition
        % let D.'*D = R.'*R where R is an upper triangular matrix
        % D.'*D should be positive definite
        R = chol(D.'*D);
        Z = R.'\(D.'*Y);
        A = R\Z;

    case 'qr'
        % Solution 3: QR decomposition
        % this method is computionally more intensive, but usually gives 
        % better numerical stablility
        [Q,R] = qr(D,0);
        A = R\(Q.'*Y);

    otherwise
        error('Unknown method! Available methods: ''inv'' (default), ''chol'' and ''qr''.')
end


Y3 = D*A; 


%%% Compute some numerical indicators of how well the fitting is
% Pearson's correlation coefficient
a = Y-mean(Y);
b = Y3-mean(Y3);
r = a.'*b / sqrt((a.'*a)*(b.'*b)); 
% root mean square error (RMSE)
Z = Y-Y3; % residuals
e = sqrt(Z.'*Z/m); 
% e2 = norm(Z); % e = e2/sqrt(m);

%%% allocate outputs
switch nargout
    case {1}
        varargout = {A};
    case {2}
        varargout = {A, Y2};
    case {3}
        varargout = {A, Y2, r};
    case {4}
        varargout = {A, Y2, r, e};
    case{0}
        % fitting (regression) result

        X_test=linspace(X(1),X(end),length(X)*100)';
        N=length(A)-1;
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

        m = length(X_test);
        D_test = cumprod([ones(m,1) X_test(:,ones(1,N))], 2) * coeff.';
        Y2 = D_test*A; 
        % no output specified, plot Y and Y2
        figure, plot(X,Y), grid on
        hold on, plot(X_test,Y2,'--r','linewidth',2)
        title(sprintf('Order N = %d; Correlation: %g; RMSE:%g', N, r, e));
        xlabel('x'), ylabel('\Sigma_{n=0}^{N} P_n(x)');
        legend('Orignial data','Regression results');
end

