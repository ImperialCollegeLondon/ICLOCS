function [const_vec_Adigator] = getAdigator4ICLOCS( X, U, t, p, data )
%genAdigator4ICLOCS - Retrieve information from Adigator for use with direct collocation method 
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
Xvec.f = X;
Xvec.dY = ones(size(Xvec.f));
Uvec.f  = U;
Uvec.dY = ones(size(Uvec.f));
Tvec.f = t;
Tvec.dY = ones(size(Tvec.f));
Pvec.f = p;
Pvec.dY = ones(size(Pvec.f));




if strcmp(data.options.ipopt.hessian_approximation,'limited-memory')
    userFunction_Y=data.data.AdigatorUserFunction_Y;
    [const_vec_Adigator] = userFunction_Y(Xvec,Uvec,Pvec,Tvec,data);
else
    userFunction_YY=data.data.AdigatorUserFunction_YY;
    [const_vec_Adigator] = userFunction_YY(Xvec,Uvec,Pvec,Tvec,data);
end


end

