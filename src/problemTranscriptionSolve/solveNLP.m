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
                    resMin_time=solution.computation_time;
%                     

                    M=data.dataNLP.sizes{7};np=data.dataNLP.sizes{2};
                    switch data.dataNLP.options.discretization
                        case{'globalLGR','hpLGR'} % p/hp Transcription Method
                            if isfield(solution,'scaledVariables')
                                [costResmin,ResColl]=costResidualMin_ModeMinRes( solution.scaledVariables.coll.X,[solution.scaledVariables.coll.U;solution.scaledVariables.coll.U(end,:)],reshape(repmat(solution.scaledVariables.coll.p,M+1,1),M+1,np),[solution.coll.T;solution.tf],data);
                            else
                                [costResmin,ResColl]=costResidualMin_ModeMinRes( solution.coll.X,[solution.coll.U;solution.coll.U(end,:)],reshape(repmat(solution.coll.p,M+1,1),M+1,np),[solution.coll.T;solution.tf],data);
                            end
                        otherwise % h Transcription Method
                            if isfield(solution,'scaledVariables')
                                [costResmin,ResColl]=costResidualMin_ModeMinRes( solution.scaledVariables.coll.X,solution.scaledVariables.coll.U,reshape(repmat(solution.scaledVariables.coll.p,M,1),M,np),solution.coll.T,data);
                            else
                                [costResmin,ResColl]=costResidualMin_ModeMinRes( solution.coll.X,solution.coll.U,reshape(repmat(solution.coll.p,M,1),M,np),solution.coll.T,data);
                            end
                    end
                    ResColl=sum(ResColl,2);
                    cost_resmin=sum(costResmin);

                    switch data.dataNLP.options.meshstrategy
                        case{'fixed','hp_flexible'}
                            ResConst=ResColl;
                            
                            idx1=ResColl<data.dataNLP.data.discErrorTol_Full;
                            ResConst(idx1)=data.dataNLP.data.discErrorTol_Full(idx1);
                            if isfield(data.dataNLP.options,'minresRelaxPct')
                                ResConst(~idx1)=(1+data.dataNLP.options.minresRelaxPct)*ResConst(~idx1);
                            else
                                ResConst(~idx1)=1.2*ResConst(~idx1);
                            end
                            if isfield(data.dataNLP.options,'minresErrorEps') 
                                ResConst(ResConst<sqrt(data.dataNLP.options.minresErrorEps))=sqrt(data.dataNLP.options.minresErrorEps);
                            end
%                             idx1=1.2*ResColl<data.dataNLP.data.discErrorTol_Full;
%                             ResConst(idx1)=ResConst(idx1)*1.2;
%                             if any(ResConst<1e-04)
%                                 ResConst(ResConst<1e-04)=ResConst(ResConst<1e-04)*10;
%                             end
%                             ResConst(ResConst<sqrt(eps))=sqrt(eps);

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
                    n=data.dataNLP.sizes{3};
                    if strcmp(data.dataNLP.options.discretization,'globalLGR') || strcmp(data.dataNLP.options.discretization,'hpLGR')
                        ng_eq=data.dataNLP.sizes{18};
                    else
                        ng_eq=data.dataNLP.sizes{15};
                    end
                    NLP.cl = [NLP.cl(1:end-n-ng_eq-1);zeros(n+ng_eq,1)];
                    NLP.cu = [NLP.cu(1:end-n-ng_eq-1);ones(n+ng_eq,1).*sqrt(ResConst)];


                    data.dataNLP.data.discErrorConstScaling=1./sqrt(ResConst)';
                    data.dataNLP.data.discErrorTol_FullScaling=data.dataNLP.data.discErrorTol_Full.*data.dataNLP.data.discErrorConstScaling';
                    ResConstScaleMat=repmat(data.dataNLP.data.discErrorConstScaling, data.nps, 1 );
                    data.ResConstScaleMat=diag(ResConstScaleMat(:));
                    data.dataNLP.options.ipopt.constr_viol_tol=min(ResConst);
                    if ~strcmp(data.dataNLP.options.NLPsolver,'NOMAD') 
                        data.dataNLP.multipliers.lambda=[solution.multipliers.lambdaNLP(1:end-n-ng_eq-1);zeros(n+ng_eq,1)];
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
                    solution.cost_resmin=cost_resmin;
                    solution.computation_time=solution.computation_time+resMin_time;
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
                    solution.cost_resmin=cost_resmin;
                    solution.computation_time=solution.computation_time+resMin_time;
                end
            elseif strcmp(data.dataNLP.options.min_res_mode,'directCostMin')

                    ResConst=data.dataNLP.data.discErrorTol_Full;
%                     data.dataNLP.data.discErrorConstScaling=1./sqrt(ResConst)';
%                     data.dataNLP.data.discErrorTol_FullScaling=data.dataNLP.data.discErrorTol_Full.*data.dataNLP.data.discErrorConstScaling';
                    ResConstScaleMat=repmat(data.dataNLP.data.discErrorConstScaling, data.nps, 1 );
                    data.ResConstScaleMat=diag(ResConstScaleMat(:));
                    data.dataNLP.options.ipopt.constr_viol_tol=min(ResConst);

                    data.mode=1;
                    data.funcs.jacobianstructure = @(data) sparse(ones(length(NLP.cl),length(NLP.z0)));
                    switch data.dataNLP.options.discretization
                        case{'globalLGR','hpLGR'} % p/hp Transcription Method
                            data.funcs.hessianstructure  = @(data) sparse(tril(ones(length(NLP.z0),length(NLP.z0))));
                        otherwise % h Transcription Method
                    end
                    [solution,status,data]=solveSingleNLP_ResidualMin(NLP,data);
                    solution.computation_time=solution.computation_time;
                
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