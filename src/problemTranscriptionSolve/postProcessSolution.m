function [solution] = postProcessSolution(z,data,solution)
%postProcessSolution - post process the solution obtained from NLP solver
%
% Syntax:  [solution] = postProcessSolution(z,data,solution)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk



if strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR')
    [nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:17});
else
    [nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:13});
end
nrc=nrcl+nrcu+nrce;


% solution.z=NLP.z0;
solution.z=z;


if strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR')
    if nt;solution.coll.tf=z(end); else solution.coll.tf=data.tf; end
    if np;solution.coll.p=z(n*(M+1)+M*m+1:n*(M+1)+M*m+np); else solution.coll.p=[]; end
    t0=z(end-nt+1);
    solution.t0=t0;

    solution.coll.X=reshape(data.map.Vx*z,M+1,n);
    solution.coll.x0=solution.coll.X(1,:)';
    solution.coll.U=reshape(data.map.Vu*z,M,m);
    if isfield(data,'t_zone_map')
        solution.coll.U=solution.coll.U(data.t_zone_map,:);
    end
    
    if data.options.adaptseg==1
        t_segment=z(end-nt+1:end)';

        tau_segment = ((1-(-1))*t_segment-(t_segment(1)+t_segment(end))) /(t_segment(end)-t_segment(1));

        tau_inc=[];

        tau_seg=cell(1,length(tau_segment)-1);
        for i=1:length(tau_segment)-1
            tau_seg{i} = ((tau_segment(i+1)+tau_segment(i)) + (tau_segment(i+1)-tau_segment(i)) * [data.map.LGR.points{npduidx(i)}; 1]) / (1-(-1));
            tau_inc=[tau_inc; tau_seg{i}(1:end-1)];
        end

        tau=diff(tau_inc);
        data.tau_inc=tau_inc;
        data.tau_segment=tau_segment;
        data.tau=tau;
    end
    
    solution.coll.T=(solution.coll.tf-t0)/2*data.tau_inc+(solution.coll.tf+t0)/2; 
    
    if ~strcmp(data.options.NLPsolver,'NOMAD') 
        solution.multipliers.lambdaNLP=solution.multipliers.lambda;
        solution.multipliers.lambdaNLP_1toN=reshape(solution.multipliers.lambda(1:M*n),M,n);
        D_N=diag(1./data.map.w);
        solution.multipliers.lambda=D_N*solution.multipliers.lambdaNLP_1toN;
        solution.multipliers.lambda_1toN=solution.multipliers.lambda;
        D_Np1=data.map.LGR.diff_matrix{end}(:,end);
        solution.multipliers.lambda(end+1,:)=D_Np1'*solution.multipliers.lambdaNLP_1toN(end-length(data.map.LGR.points{end})+1:end,:);
        
        if ng
            solution.multipliers.lambda_g=zeros(ng,M);
            solution.multipliers.lambda_g(logical(data.gActiveIdx'))=solution.multipliers.lambdaNLP(n*M+1:n*M+ngActive);
            solution.multipliers.lambda_g=solution.multipliers.lambda_g';
        end
        if nrc
            solution.multipliers.lambda_rc=solution.multipliers.lambdaNLP((M*(n)+ngActive+1):M*(n)+ngActive+nrc);
        end
        if nb
            solution.multipliers.lambda_b=solution.multipliers.lambdaNLP((M*(n)+ngActive+nrc+1):end);
        end
    end
    
        solution.X=solution.coll.X;
        solution.U=solution.coll.U;
        solution.T=solution.coll.T;
        solution.p=solution.coll.p;
        solution.x0=solution.coll.x0;
        solution.tf=solution.coll.tf;
else
    
    if nt==1
        solution.coll.t0=data.t0;
        solution.coll.tf=z(1);
        t0=data.t0;
    elseif nt==2
        t0=z(1);
        tf=z(2);
        solution.coll.t0=t0;    
        solution.coll.tf=tf;
    else
        t0=data.t0;
        solution.coll.t0=data.t0;   
        solution.coll.tf=data.tf; 
    end
    
    if np
        solution.coll.p=z(nt+1:nt+np);
    else
        solution.coll.p=[]; 
    end


    if strcmp(data.options.transcription,'multiple_shooting')


        solution.X=reshape(data.map.Vx*z,n,M)';
        solution.x0=solution.X(1,:)';
        usp=data.map.Vu*z;
        if isfield(data,'t_zone_map')
            usp=usp(data.t_zone_map,:);
        end
        solution.U=reshape([usp;usp(end-m+1:end)],m,M)';
        solution.T=(solution.coll.tf-t0)*[0;data.tau_inc]+t0; 

    else

        solution.coll.X=reshape(data.map.Vx*z,n,M)';
        solution.coll.x0=solution.coll.X(1,:)';
        usp=reshape(data.map.Vu*z,m,N)';
        if isfield(data,'t_zone_map')
            usp=usp(data.t_zone_map,:);
        end
        solution.coll.U=kron(usp,ones((M)/N,1));
        solution.coll.T=(solution.coll.tf-solution.coll.t0)*[0;data.tau_inc*data.Nm/ns]+solution.coll.t0;

        if ~strcmp(data.options.NLPsolver,'NOMAD') 
            solution.multipliers.lambdaNLP=solution.multipliers.lambda;
            solution.multipliers.lambdaNLP_1toM=reshape(solution.multipliers.lambdaNLP(1:n*M)',n,M)';
        end
        
        if strcmp(data.options.discretization,'hermite')
            solution.X=solution.coll.X(1:2:end,:);
            solution.U=solution.coll.U(1:2:end,:);
            solution.T=solution.coll.T(1:2:end);
            if ~strcmp(data.options.NLPsolver,'NOMAD') 
                solution.multipliers.lambdaNLP_1toN=solution.multipliers.lambdaNLP_1toM(1:2:end,:);
                solution.multipliers.lambda_1toN=solution.multipliers.lambdaNLP_1toN;
            end
        else
            solution.X=solution.coll.X;
            solution.U=solution.coll.U;
            solution.T=solution.coll.T;
            solution.multipliers.lambdaNLP_1toN=solution.multipliers.lambdaNLP_1toM;
        end
        solution.p=solution.coll.p;
        solution.x0=solution.coll.x0;
        solution.t0=solution.coll.t0;
        solution.tf=solution.coll.tf;

        if ~strcmp(data.options.NLPsolver,'NOMAD') 
            if ng
                solution.multipliers.lambda_g=zeros(ng,M);
                solution.multipliers.lambda_g(logical(data.gActiveIdx'))=solution.multipliers.lambdaNLP(n*M+1:n*M+ngActive);
                solution.multipliers.lambda_g=solution.multipliers.lambda_g';
            end
            if nrc
                solution.multipliers.lambda_rc=solution.multipliers.lambdaNLP(n*M+ngActive+1:n*M+ngActive+nrc);
            end
            if nb
                solution.multipliers.lambda_b=solution.multipliers.lambdaNLP(n*M+ngActive+nrc+1:n*M+ngActive+nrc+nb);
            end
        end
    end
end

%scaling back
if data.options.scaling
    solution.scaledVariables.X=solution.X;
    solution.scaledVariables.x0=solution.x0;
    solution.scaledVariables.U=solution.U;
    solution.scaledVariables.coll.X=solution.coll.X;
    solution.scaledVariables.coll.x0=solution.coll.x0;
    solution.scaledVariables.coll.U=solution.coll.U;
    solution.X=scale_variables_back( solution.X, data.data.Xscale_back, data.data.Xshift );
    solution.x0=scale_variables_back( solution.x0', data.data.Xscale_back, data.data.Xshift )';
    solution.U=scale_variables_back( solution.U, data.data.Uscale_back, data.data.Ushift );
    solution.coll.X=scale_variables_back( solution.coll.X, data.data.Xscale_back, data.data.Xshift );
    solution.coll.x0=scale_variables_back( solution.coll.x0', data.data.Xscale_back, data.data.Xshift )';
    solution.coll.U=scale_variables_back( solution.coll.U, data.data.Uscale_back, data.data.Ushift );
    if isfield(data.data,'Pscale')
        solution.scaledVariables.p=solution.p;
        solution.scaledVariables.coll.p=solution.coll.p;
        solution.p=scale_variables_back( solution.p', data.data.Pscale_back, data.data.Pshift )';
        solution.coll.p=scale_variables_back( solution.coll.p', data.data.Pscale_back, data.data.Pshift )';
    else
        solution.scaledVariables.p=[];
        solution.scaledVariables.coll.p=[];
    end
end

end

