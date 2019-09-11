function [ ConstraintError ] = calcConstraintViolation( data, guess, x_guess, u_guess, problem, T)
%calcConstraintViolation - check the extend of contraint violation for path
%constraints in the provided initial guess
%
% Syntax:   [ ConstraintError ] = calcConstraintViolation( data, guess, x_guess, u_guess, problem, T)
%
% Inputs:
%    data, guess, problem - Defined in transcribeOCP
%    x_guess, u_guess - the inital guesses
%    T - the time vector
%
% Outputs:
%    ConstraintError - the contraint violation error for the initial guess
%     
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

g=data.functions_unscaled{4};

G=g(x_guess,u_guess,guess.parameters,T,data.data);           % Function evaluations.
Gl=G(:,~problem.data.gFilter)-problem.constraints.gl;
Gl(Gl>0)=0;
Gu=problem.constraints.gu-G(:,~problem.data.gFilter);
Gu(Gu>0)=0;
    
ConstraintError=[abs(Gu),abs(Gl)];
end

