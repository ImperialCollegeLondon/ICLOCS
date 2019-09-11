function [solution,status,data]=solveNLP(NLP,data)
%SOLVENLP - Solve the nonlinear program and return the solution
%
% Syntax:  solution = solveNLP(NLP,options)
%
% Inputs:
%    NLP     - Information required by the NLP solver
%    data    - Data passed to the functions evaluated during optimization
%
% Outputs:
%    solution - Structure containing the optimal final time, parameters,
%               states, controls and Lagrange multipliers.
%    status - variable containing the exit condition of the solver. See the
%             fmincon and ipopt description 
%
% Other m-files required: ipopt.mex and or fmincon.m
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

% Solve the nonlinear program:

global mode_min_res


if isfield(data,'mpdata')
    if ~isfield(data,'dataNLP')
        if isfield(data.mpdata.data,'penalty') && strcmp(data.mpdata.options.mp.regstrategy,'reg_priority')
            if isfield(data.mpdata.data.penalty,'i') && isfield(data.mpdata.data.penalty,'values')
                for i=1:length(data.mpdata.data.penalty.values)
                    data.mpdata.data.penalty.i=i;
                    for j=1:length(data.phasedata)
                        data.phasedata{j}.data.penalty.i=i;
                    end
                    [solution,status,data]=solveSingleNLP_DirectCollocation_MultiPhase(NLP,data);
                    for j=1:length(solution.phaseSol)
                        data.phasedata{j}.multipliers.lambda=solution.phaseSol{j}.multipliers.lambdaNLP;
                    end
                    NLP.mpinfoNLP.z0=solution.mp.z_org;
                end
            else
                error('Regularization Parameters Not Properly Configured!')
            end
        else
            [solution,status,data]=solveSingleNLP_DirectCollocation_MultiPhase(NLP,data);
        end
        
    end
    
    
else
    

    if isfield(data,'dataNLP')

        if isfield(data.dataNLP.data,'penalty') && strcmp(data.dataNLP.options.regstrategy,'reg_priority')
            if isfield(data.dataNLP.data.penalty,'i') && isfield(data.dataNLP.data.penalty,'values')
                for i=1:length(data.dataNLP.data.penalty.values)
                    data.dataNLP.data.penalty.i=i;
                    [solution,status,data]=solveSingleNLP_ResidualMin(NLP,data);
                    data.dataNLP.multipliers.lambda=solution.multipliers.lambdaNLP;
                    NLP.z0=solution.z;
                end
            else
                error('Regularization Parameters Not Properly Configured!')
            end
        else
            if strcmp(data.dataNLP.options.min_res_mode,'alternating')
                    data.mode=2;
                    if data.dataNLP.options.resminEarlyStop
                        data.funcs.iterfunc=@callback_minres;
                    end
                    [solution,status,data]=solveSingleNLP_ResidualMin(NLP,data);
                    if data.dataNLP.options.resminEarlyStop
                        data.funcs=rmfield(data.funcs,'iterfunc');
                    end

                    M=data.dataNLP.sizes{7};np=data.dataNLP.sizes{2};
                    switch data.dataNLP.options.discretization
                        case{'globalLGR','hpLGR'} % p/hp Transcription Method
                            if isfield(solution,'scaledVariables')
                                [~,ResColl]=costResidualMin_ModeMinRes( solution.scaledVariables.coll.X,[solution.scaledVariables.coll.U;solution.scaledVariables.coll.U(end,:)],reshape(repmat(solution.scaledVariables.coll.p,M+1,1),M+1,np),[solution.coll.T;solution.tf],data);
                            else
                                [~,ResColl]=costResidualMin_ModeMinRes( solution.coll.X,[solution.coll.U;solution.coll.U(end,:)],reshape(repmat(solution.coll.p,M+1,1),M+1,np),[solution.coll.T;solution.tf],data);
                            end
                            ResColl=sum(ResColl,2);
                        otherwise % h Transcription Method
                            if isfield(solution,'scaledVariables')
                                [~,ResColl]=costResidualMin_ModeMinRes( solution.scaledVariables.coll.X,solution.scaledVariables.coll.U,reshape(repmat(solution.scaledVariables.coll.p,M,1),M,np),solution.coll.T,data);
                            else
                                [~,ResColl]=costResidualMin_ModeMinRes( solution.coll.X,solution.coll.U,reshape(repmat(solution.coll.p,M,1),M,np),solution.coll.T,data);
                            end
                            ResColl=sum(ResColl,2);
                    end

                    switch data.dataNLP.options.meshstrategy
                        case{'fixed','hp_flexible'}
                            ResConst=ResColl;
                            ResConst=ResConst*2;
                            ResConst(ResConst<sqrt(eps))=sqrt(eps);

                        case{'mesh_refinement'}
                            if data.dataNLP.options.scaling
                                ResConst=ResColl;
                                ResConst(ResConst<sqrt(eps))=sqrt(eps);

                            else
                                switch data.dataNLP.options.errortype
                                    case{'local_abs'}
                                        ResConst=ResColl;
                                        ResConst(ResConst<sqrt(eps))=sqrt(eps);
                                    case{'int_res'}
                                        ResConst=max(ResColl,data.dataNLP.data.discErrorTol_Full/2);
                                        ResConst(ResConst<sqrt(eps))=sqrt(eps);
                                    case{'both'}
                                        ResConst=min(ResColl*2,max(ResColl,data.dataNLP.data.discErrorTol_Full));
                                        ResConst(ResConst<sqrt(eps))=sqrt(eps);
                                end
                            end
                    end

                    NLP.z0 = solution.z;
                    NLP.cl = [NLP.cl(1:end-data.dataNLP.sizes{3}-1);zeros(data.dataNLP.sizes{3},1)];
                    NLP.cu = [NLP.cu(1:end-data.dataNLP.sizes{3}-1);ones(data.dataNLP.sizes{3},1)];


                    data.dataNLP.data.discErrorConstScaling=1./ResConst';
                    data.dataNLP.data.discErrorTol_FullScaling=data.dataNLP.data.discErrorTol_Full.*data.dataNLP.data.discErrorConstScaling';
                    ResConstScaleMat=repmat(data.dataNLP.data.discErrorConstScaling, data.nps, 1 );
                    data.ResConstScaleMat=diag(ResConstScaleMat(:));
                    data.dataNLP.options.ipopt.constr_viol_tol=min(ResConst);
                    if ~strcmp(data.dataNLP.options.NLPsolver,'NOMAD') 
                        data.dataNLP.multipliers.lambda=[solution.multipliers.lambdaNLP(1:end-data.dataNLP.sizes{3}-1);zeros(data.dataNLP.sizes{3},1)];
                    end

                if (all(ResColl<data.dataNLP.data.discErrorTol_Full) && strcmp(data.dataNLP.options.errortype,'int_res')) || strcmp(data.dataNLP.options.meshstrategy,'fixed')
                    data.mode=1;
                    data.funcs.jacobianstructure = @(data) sparse(ones(length(NLP.cl),length(NLP.z0)));
                    switch data.dataNLP.options.discretization
                        case{'globalLGR','hpLGR'} % p/hp Transcription Method
                            data.funcs.hessianstructure  = @(data) sparse(tril(ones(length(NLP.z0),length(NLP.z0))));
                        otherwise % h Transcription Method
                    end
                    [solution,status,data]=solveSingleNLP_ResidualMin(NLP,data);
                else
                    if mode_min_res && ~strcmp(data.dataNLP.options.errortype,'int_res')
                        data.mode=1;
                        data.funcs.jacobianstructure = @(data) sparse(ones(length(NLP.cl),length(NLP.z0)));
                        switch data.dataNLP.options.discretization
                            case{'globalLGR','hpLGR'} % p/hp Transcription Method
                                data.funcs.hessianstructure  = @(data) sparse(tril(ones(length(NLP.z0),length(NLP.z0))));
                            otherwise % h Transcription Method
                        end
                        [solution,status,data]=solveSingleNLP_ResidualMin(NLP,data);
                    end

                    if mode_min_res && all(ResColl<data.dataNLP.data.discErrorTol_Full)
                        mode_min_res=0;
                        solution.min_res_satisfy=1;
                    end
                end
            elseif strcmp(data.dataNLP.options.min_res_mode,'weightedCost')
                data.mode=0;
                [solution,status,data]=solveSingleNLP_ResidualMin(NLP,data);
            else
                error('Integral Residual Minimization Optimization Method Not Properly Chosen!')       
            end
        end
    else

        if isfield(data.data,'penalty') && strcmp(data.options.regstrategy,'reg_priority')
            if isfield(data.data.penalty,'i') && isfield(data.data.penalty,'values')
                for i=1:length(data.data.penalty.values)
                    data.data.penalty.i=i;
                    [solution,status,data]=solveSingleNLP_DirectCollocation(NLP,data);
                    data.multipliers.lambda=solution.multipliers.lambdaNLP;
                    NLP.z0=solution.z;
                end
            else
                error('Regularization Parameters Not Properly Configured!')
            end
        else
            [solution,status,data]=solveSingleNLP_DirectCollocation(NLP,data);
        end

    end

end