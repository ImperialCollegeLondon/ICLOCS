function [options,guess] = prepareReStart(options,guess,solution,OCP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isfield(options,'mp')
    OCP.infoNLP.mpinfoNLP.z0=solution.mp.z_org;
    OCP.data.mpdata.multipliers.lambda=solution.mp.multipliers.lambda;
else
    [ options, guess] = doWarmStart( options, guess, solution, OCP.data );
end

