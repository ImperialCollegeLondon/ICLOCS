function [options,guess] = prepareWarmStart(options,guess,solution,OCP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isfield(options,'mp')
    nphase=length(solution.phaseSol);
    for j=1:nphase
        [ options.phaseoptions{j}, guess.phases{j}] = doWarmStart( OCP.options.phaseoptions{j}, guess.phases{j}, solution.phaseSol{j}, OCP.data.phasedata{j} );
        guess.mp.lambda_nbl=solution.mp.multipliers.lambda(end-OCP.data.mpdata.mpsizes.nbl_l+1:end);
    end
else
    [ options, guess] = doWarmStart( options, guess, solution, OCP.data );
end

