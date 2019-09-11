function [infoNLP,data,options]=transcribeOCP(problem,guess,options)
%TRANSCRIBEOCP - Process information from 'problem', 'guess' and 'options' for NLP solver
%Specifically:
%Error checking of function definitions and bounds
%Define bounds for NLP variable + continuity, path and boundary constraints
%Format matrices for direct transcription method(if required)
%Generate initial guess for optimization
%Generate structure of the jacobian of the constraints(if required)
%Construct optimal finite-difference pertubation sets(if required)
%
% Syntax:  [infoNLP,data]=transcribeOCP(problem,guess,options)
%
% Inputs:
%    problem - Optimal control problem definition
%    guess   - Guess used to generate starting point for optimization
%    options - Settings in file settings.m
%
% Outputs:
%    infoNLP - Information required by the NLP solver
%    data - Data passed to the functions evaluated during optimization
%    options - Settings after processing
%
% Other m-files required: scale_problem, getStructure, getStructureA, getPertubations.
% Subfunctions: checkErrors, transcriptionMatrix
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

global sol; 
sol=[];                             % Initialize solution structure

% Running mode
if ~isfield(problem.data,'mode')
    problem.data.mode.currentMode='Original';
end

[infoNLP,data,options]=transcribeOCP_eachPhase(problem,guess,options);
