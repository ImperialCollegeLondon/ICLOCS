function [ ResNorm_vec,Res_vec ] = getAdigator4ICLOCS_minres( X, U, T, P, data )
%genAdigator4ICLOCS - Retrieve information from Adigator for use with integrated residual minimization method
%
% Syntax:  [const_vec_Adigator] = getAdigator4ICLOCS( X, U, t, p, data )
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

dataNLP=data.dataNLP;
switch dataNLP.options.discretization
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        Xvec.f = X;
        Xvec.dY = ones(size(Xvec.f));
        Uvec.f  = U(1:end-1,:);
        Uvec.dY = ones(size(Uvec.f));
        Tvec.f = T;
        Tvec.dY = ones(size(Tvec.f));
        Pvec.f = P;
        Pvec.dY = ones(size(Pvec.f));
    otherwise % h Transcription Method
        Xvec.f = X;
        Xvec.dY = ones(size(Xvec.f));
        Uvec.f  = U;
        Uvec.dY = ones(size(Uvec.f));
        Tvec.f = T;
        Tvec.dY = ones(size(Tvec.f));
        Pvec.f = P;
        Pvec.dY = ones(size(Pvec.f));
end

if strcmp(dataNLP.options.ipopt.hessian_approximation,'limited-memory')
    [ResNorm_vec,Res_vec] = minresCost_Y(Xvec,Uvec,Pvec,Tvec,data);
else
    [ResNorm_vec,Res_vec] = minresCost_YY(Xvec,Uvec,Pvec,Tvec,data);
end


end

