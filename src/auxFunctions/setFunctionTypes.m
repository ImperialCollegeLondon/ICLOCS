function [ data ] = setFunctionTypes( problem, data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isfield(problem,'FcnTypes') 
    if (strcmp(problem.FcnTypes.StageCost,'Constant')) || (strcmp(problem.FcnTypes.StageCost,'None'))
        data.FD.FcnTypes.Ltype=0;
    elseif (strcmp(problem.FcnTypes.StageCost,'Linear'))
        data.FD.FcnTypes.Ltype=2;
    else
        data.FD.FcnTypes.Ltype=1;
    end
    
    if (strcmp(problem.FcnTypes.TerminalCost,'Constant')) || (strcmp(problem.FcnTypes.TerminalCost,'None'))
        data.FD.FcnTypes.Etype=0;
    elseif (strcmp(problem.FcnTypes.TerminalCost,'Linear'))
        data.FD.FcnTypes.Etype=2;
    else
        data.FD.FcnTypes.Etype=1;
    end
    
    if (strcmp(problem.FcnTypes.TerminalConst,'Constant')) || (strcmp(problem.FcnTypes.TerminalConst,'None'))
        data.FD.FcnTypes.Btype=0;
    elseif (strcmp(problem.FcnTypes.TerminalConst,'Linear'))
        data.FD.FcnTypes.Btype=2;
    else
        data.FD.FcnTypes.Btype=1;
    end
else
    data.FD.FcnTypes.Ltype=1;
    data.FD.FcnTypes.Etype=1;
    data.FD.FcnTypes.Btype=1;
end

% problem.FcnTypes.Dynamics='Nonlinear';
% problem.FcnTypes.PathConstraint='None';
% problem.FcnTypes.StageCost='Constant';
% problem.FcnTypes.TerminalCost='Constant';
% problem.FcnTypes.TerminalConst='None';
end

