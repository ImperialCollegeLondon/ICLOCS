function [ options, guess] = doMeshRefinement( options, problem, guess, data, solution, mriteration )
%doMeshRefinement - Automatically run meshrefinement scheme
%
% Syntax:  [ options, guess] = doMeshRefinement( options, problem, guess, data, solution, mriteration )
%
% Inputs:
%    options, problem, guess, data - Defined in transcribeOCP
%    solutions - Structure containing the solution
%    minItervalScale - mesh refinement iteration number
%
% Outputs:
%    options, guess - files for next iteration of OCP run, with warm
%    starting
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


if isfield(data,'dataNLP')
    data=data.dataNLP;
end

options_old=options;

minItervalScale=0.9^(mriteration-1);
xu=problem.states.xu;xu(isinf(xu))=500;
xl=problem.states.xl;xl(isinf(xl))=-500;
Xscales=1./(xu-xl);

discErrorTol=options.discErrorTol;
discErrorTol_Scaled= scale_variables( discErrorTol, Xscales, 0 );
options.discErrorTol_Scaled= discErrorTol_Scaled;
solution.ErrorMax = max( solution.Error, [], 1 );
solution.ErrorScaled = scale_variables( solution.Error, Xscales, 0 );
error_ratio=solution.ErrorMax./discErrorTol;
nt=data.sizes{1};

if strcmp(options.MRstrategy,'efficient')
    MeshRefinement_h=@MeshRefinement_h_OL;
    MeshRefinement_Auto=@MeshRefinement_Auto_OL;
elseif strcmp(options.MRstrategy,'aggressive')
    MeshRefinement_h=@MeshRefinement_h_TO;
    MeshRefinement_Auto=@MeshRefinement_Auto_TO;
else
    disp('selected mesh refinement priority not recognized, will default to aggressive')
    MeshRefinement_h=@MeshRefinement_h_TO;
    MeshRefinement_Auto=@MeshRefinement_Auto_TO;
end

settings=problem.settings;

if isfield(options,'AutoDirect')
    const_active_ratio=solution.NumActiveConstraint/data.sizes{3};
    if strcmp(options.discretization,'hpLGR') || strcmp(options.discretization,'globalLGR')
        fcn_smoothness=max(std(abs([diff(solution.X(1:end-1,:))./repmat(diff(solution.T),1,size(solution.X,2)) diff(solution.U)./repmat(diff(solution.T),1,size(solution.U,2))]))./median(abs([diff(solution.X(1:end-1,:))./repmat(diff(solution.T),1,size(solution.X,2)) diff(solution.U)./repmat(diff(solution.T),1,size(solution.U,2))])));
    else
        fcn_smoothness=max(std(abs([diff(solution.X)./repmat(diff(solution.T),1,size(solution.X,2)) diff(solution.U)./repmat(diff(solution.T),1,size(solution.U,2))]))./median(abs([diff(solution.X)./repmat(diff(solution.T),1,size(solution.X,2)) diff(solution.U)./repmat(diff(solution.T),1,size(solution.U,2))])));
    end
    if fcn_smoothness<=10 && max(error_ratio)<=10 && ~isfield(problem.inputs,'url')
        options= settings(1,40);
        options.constraintErrorTol_org=problem.constraintErrorTol;
    elseif const_active_ratio<=0.2 && ~isfield(problem.inputs,'url')
        options= settings(1,8);
        options.constraintErrorTol_org=problem.constraintErrorTol;
    else
        options= settings(40);
        options.discretization='hermite';
        options.resultRep=options_old.resultRep;
        options.discErrorTol=discErrorTol;
        options.constraintErrorTol=problem.constraintErrorTol;
        options.discErrorTol_Scaled= discErrorTol_Scaled;
        options.constraintErrorTol_org=problem.constraintErrorTol;
        if options.tau==0
            M=deal(data.sizes{7});
            ns=deal(data.sizes{9});
            options.tau=diff(linspace(0,1,M/ns+1));
            [ options ] = MeshRefinement_h( options, data, solution ); 
        else
            [ options ] = MeshRefinement_h( options, data, solution );
        end
    end
    options.discErrorTol=discErrorTol;
    options.constraintErrorTol=problem.constraintErrorTol;
    options.discErrorTol_Scaled= discErrorTol_Scaled;
else
    if isfield(options,'hpAdaptive')
        t_seg=solution.z((end-nt+1):end)';
        tau_seg=2*(t_seg-min(t_seg))/(max(t_seg)-min(t_seg))-1;
        options= settings(8*ones(1,length(tau_seg)-1),tau_seg);
        options.constraintErrorTol_org=problem.constraintErrorTol;
        options.resultRep='default';
        options.discErrorTol=discErrorTol;
        options.constraintErrorTol=problem.constraintErrorTol;
        options.discErrorTol_Scaled= discErrorTol_Scaled;
    end
    if strcmp(options.discretization,'hpLGR')
        if (strcmp(options.MeshRefinement,'Auto')) 
            [ options ] = MeshRefinement_Auto( options, data, solution, minItervalScale );
        elseif (strcmp(options.MeshRefinement,'AI')) 
            error('For current version, use Auto for mesh refinement')
            [ options ] = MeshRefinement_AI( options, data, solution, minItervalScale );
        elseif (strcmp(options.MeshRefinement,'IO')) 
            error('For current version, use Auto for mesh refinement')
            [ options ] = MeshRefinement_IO( options, data, solution );
        end
    elseif strcmp(options.discretization,'globalLGR')
        options.pdegree=options.pdegree+2;
        if options.pdegree>=40
            if max(error_ratio)>=50 
                options= settings(40);
                options.constraintErrorTol_org=problem.constraintErrorTol;
                options.resultRep='default';
            else
                options= settings(5,8);
                options.constraintErrorTol_org=problem.constraintErrorTol;
                options.resultRep='default';
            end
            options.discErrorTol=discErrorTol;
            options.constraintErrorTol=problem.constraintErrorTol;
            options.discErrorTol_Scaled= discErrorTol_Scaled;
        end
    else
        if options.tau==0
            M=deal(data.sizes{7});
            ns=deal(data.sizes{9});
            options.tau=diff(linspace(0,1,M/ns+1));
            [ options ] = MeshRefinement_h( options, data, solution ); 
        else
            [ options ] = MeshRefinement_h( options, data, solution );
        end
    end
end

[ options, guess] = doWarmStart( options, guess, solution, data );
    




