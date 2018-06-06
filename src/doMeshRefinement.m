function [ options, guess] = doMeshRefinement( options, problem, guess, data, solution, mriteration )
%doMeshRefinement - Automatically run meshrefinement scheme
%
% Syntax:  [ options, guess] = doMeshRefinement( options, problem, guess, data, solution, minItervalScale)
%
% Inputs:
%    options, problem, guess, data - Defined in transcribeOCP
%    solutions - Structure containing the solution
%    minItervalScale - scaling of the minimum interval length
%
% Outputs:
%    options, guess - files for next iteration of OCP run, with warm
%    starting
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

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

if isfield(options,'AutoDirect')
    const_active_ratio=solution.NumActiveConstraint/data.sizes{3};
    if strcmp(options.transcription,'hpLGR') || strcmp(options.transcription,'globalLGR')
        fcn_smoothness=max(std(abs([diff(solution.X(1:end-1,:))./repmat(diff(solution.T),1,size(solution.X,2)) diff(solution.U)./repmat(diff(solution.T),1,size(solution.U,2))]))./median(abs([diff(solution.X(1:end-1,:))./repmat(diff(solution.T),1,size(solution.X,2)) diff(solution.U)./repmat(diff(solution.T),1,size(solution.U,2))])));
    else
        fcn_smoothness=max(std(abs([diff(solution.X)./repmat(diff(solution.T),1,size(solution.X,2)) diff(solution.U)./repmat(diff(solution.T),1,size(solution.U,2))]))./median(abs([diff(solution.X)./repmat(diff(solution.T),1,size(solution.X,2)) diff(solution.U)./repmat(diff(solution.T),1,size(solution.U,2))])));
    end
    if fcn_smoothness<=10 && max(error_ratio)<=10 && ~isfield(problem.inputs,'url')
        options= settings_hp(1,40);
        options.resultRep='default';
    elseif const_active_ratio<=0.2 && fcn_smoothness<=10 && ~isfield(problem.inputs,'url')
        options= settings_hp(1,8);
        options.resultRep='default';
    else
        options= settings_h(40);
        options.resultRep='default';
        options.discErrorTol=discErrorTol;
        options.constraintErrorTol=problem.constraintErrorTol;
        options.discErrorTol_Scaled= discErrorTol_Scaled;
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
        options= settings_hp(8*ones(1,length(tau_seg)-1),tau_seg);
        options.resultRep='default';
        options.discErrorTol=discErrorTol;
        options.constraintErrorTol=problem.constraintErrorTol;
        options.discErrorTol_Scaled= discErrorTol_Scaled;
    end
    if strcmp(options.transcription,'hpLGR')
        if (strcmp(options.MeshRefinement,'Auto')) 
            [ options ] = MeshRefinement_Auto( options, data, solution, minItervalScale );
        elseif (strcmp(options.MeshRefinement,'AI')) 
            [ options ] = MeshRefinement_AI( options, data, solution, minItervalScale );
        elseif (strcmp(options.MeshRefinement,'IO')) 
            [ options ] = MeshRefinement_IO( options, data, solution );
        end
    elseif strcmp(options.transcription,'globalLGR')
        options.pdegree=options.pdegree+2;
        if options.pdegree>=40
            if max(error_ratio)>=50 
                options= settings_h(40);
                options.resultRep='default';
            else
                options= settings_hp(5,8);
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

% if strcmp(data.options.transcription,'multiple_shooting')
%     options.start='Cold';
% else
    if (isfield(options,'hpAdaptive') && options.hpAdaptive==1)
        options.start='Warm';
    else
        options.start='Hot';
    end
    guess.states=solution.X;
    guess.inputs=solution.U;
    guess.time=solution.T;
    guess.tf=solution.tf;
    guess.multipliers.lambda=solution.multipliers.lambda_1toN;
    ng=deal(data.sizes{5});
    nb=deal(data.sizes{6});
    if ng
        if strcmp(data.options.transcription,'hermite') || strcmp(data.options.transcription,'AutoDirect')
            guess.timeFull=solution.org.T;
            guess.multipliers.lambda_g=solution.multipliers.lambda_g;
        else
            guess.timeFull=solution.T;
            guess.multipliers.lambda_g=solution.multipliers.lambda_g;
        end
    end
    if nb
        guess.multipliers.lambda_b=solution.multipliers.lambda_b;
    end

% end




