function [grad]=resMin_gradientCost_ModeMinRes(X,U,P,T,t0,tf,data)
%resMin_gradientCost_ModeMinRes - cost gradient computation for
%integrated residual minimization (alternating method: residual
%minimization) using finite difference
%
% Syntax:   [grad]=resMin_gradientCost_ModeMinRes(X,U,P,T,t0,tf,data)
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

dataNLP=data.dataNLP;

e=dataNLP.options.perturbation.J;                                 % Pertubation size
[nt,np,n,m,ng,nb,M,N]=deal(dataNLP.sizes{1:8});

% Compute contribution to the gradient
idx=dataNLP.FD.index.Ly;
nfd=size(idx,2);                               
et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;


switch dataNLP.options.discretization
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        nx=(M+1)*n;                               % Number of unknown states
        nu=M*m;                               % Number of unknown controls

        
        if data.free_time
            i_st=1;
            i_end=nfd;
            nz=nx+nu+nt+np;                       % Length of the optimization variable 
        else
            i_st=1;
            i_end=m+n+np;
            nz=nx+nu+np;                       % Length of the optimization variable 
        end
        
        Resz=zeros(1,nz);
        tau_seg_idx=[dataNLP.tau_seg_idx;dataNLP.tau_seg_idx(end)];
        idx_state=[idx(:,1:n);idx(end,1:n)+1];
        for i=i_st:i_end
            if (data.free_time && i<=(m+n)) || (~data.free_time)
                for j=1:size(data.idx_perturb,2)
                    if i<=n
                        [ResNormCost_p,~]=costResidualMin_ModeMinRes( X+[ex{i};ex{i}(end,:)].*data.idx_perturb(:,j)*e,U+[eu{i};eu{i}(end,:)].*data.idx_perturb(:,j)*e,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*[T;1]+(tf+etf{i}*e+t0+et0{i}*e)/2,data);
                        [ResNormCost_m,~]=costResidualMin_ModeMinRes( X-[ex{i};ex{i}(end,:)].*data.idx_perturb(:,j)*e,U-[eu{i};eu{i}(end,:)].*data.idx_perturb(:,j)*e,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*[T;1]+(tf-etf{i}*e+t0-et0{i}*e)/2,data);
                        dRes=(ResNormCost_p-ResNormCost_m)/(2*e);
                    else
                        [ResNormCost_p,~]=costResidualMin_ModeMinRes( X+[ex{i};ex{i}(end,:)].*data.idx_perturb(:,j)*e,U+[eu{i};eu{i}(end,:)].*[data.idx_perturb(1:end-1,j);data.idx_perturb(end-1,j)]*e,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*[T;1]+(tf+etf{i}*e+t0+et0{i}*e)/2,data);
                        [ResNormCost_m,~]=costResidualMin_ModeMinRes( X-[ex{i};ex{i}(end,:)].*data.idx_perturb(:,j)*e,U-[eu{i};eu{i}(end,:)].*[data.idx_perturb(1:end-1,j);data.idx_perturb(end-1,j)]*e,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*[T;1]+(tf-etf{i}*e+t0-et0{i}*e)/2,data);
                        dRes=(ResNormCost_p-ResNormCost_m)/(2*e);
                    end
                    if j==1
                         if mod(data.nps,2)
                            dRes=[dRes(1) dRes(2:2:end)+dRes(3:2:end)];
                         else
                            dRes=[dRes(1) dRes(2:2:end)+[dRes(3:2:end) 0]];
                         end
                         if i<=n
                            Resz=Resz+sparse(1,idx_state(logical(data.idx_perturb(:,j)),i),dRes,1,nz);
                         else
                            Resz=Resz+sparse(1,idx(logical(data.idx_perturb(1:end-1,j)),i),dRes(1:length(idx(logical(data.idx_perturb(1:end-1,j)),i))),1,nz);
                         end
                    elseif j==2
                         if mod(data.nps,2)
                            dRes=dRes(1:2:end)+[dRes(2:2:end) 0];
                         else
                            dRes=dRes(1:2:end)+dRes(2:2:end);
                         end
                         if i<=n
                            Resz=Resz+sparse(1,idx_state(logical(data.idx_perturb(:,j)),i),dRes,1,nz);
                         else
                            Resz=Resz+sparse(1,idx(logical(data.idx_perturb(1:end-1,j)),i),dRes(1:length(idx(logical(data.idx_perturb(1:end-1,j)),i))),1,nz);
                         end
                    else
                         dRes=dRes(tau_seg_idx(logical(data.idx_perturb(:,j))));
                         if i<=n
                            Resz=Resz+sparse(1,idx_state(logical(data.idx_perturb(:,j)),i),dRes,1,nz);
                         else
                            Resz=Resz+sparse(1,idx(logical(data.idx_perturb(1:end-1,j)),i),dRes(1:length(idx(logical(data.idx_perturb(1:end-1,j)),i))),1,nz);
                         end
                    end
                end
            else
                [ResNormCost_p,~]=costResidualMin_ModeMinRes( X,U,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*[T;1]+(tf+etf{i}*e+t0+et0{i}*e)/2,data);
                [ResNormCost_m,~]=costResidualMin_ModeMinRes( X,U,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*[T;1]+(tf-etf{i}*e+t0-et0{i}*e)/2,data);
                dRes=sum((ResNormCost_p-ResNormCost_m)/(2*e));
                Resz=Resz+sparse(1,idx(1,i),dRes,1,nz);
            end
        end
    otherwise % h Transcription Method
        nx=M*n;                               % Number of unknown states
        nu=N*m;                               % Number of unknown controls

        if nt
            i_st=1;
            i_end=nfd;
            nz=nx+nu+nt+np;                       % Length of the optimization variable 
        else
            i_st=nt+1;
            i_end=nt+np+m+n;
            idx=idx-nt;
            nz=np+nx+nu;                       % Length of the optimization variable 
        end
        Resz=zeros(1,nz);
        
         for i=i_st:i_end
            if (nt && i>(nt+np)) || (~nt)
                for j=1:size(data.idx_perturb,2)
                    [ResNormCost_p,~]=costResidualMin_ModeMinRes( X+ex{i}.*data.idx_perturb(:,j)*e,U+eu{i}.*data.idx_perturb(:,j)*e,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)*T+t0+et0{i}*e,data);
                    [ResNormCost_m,~]=costResidualMin_ModeMinRes( X-ex{i}.*data.idx_perturb(:,j)*e,U-eu{i}.*data.idx_perturb(:,j)*e,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)*T+t0-et0{i}*e,data);
                    dRes=(ResNormCost_p-ResNormCost_m)/(2*e);
                    if j==1
                         if mod(data.nps,2)
                            dRes=[dRes(1) dRes(2:2:end)+dRes(3:2:end)];
                         else
                            dRes=[dRes(1) dRes(2:2:end)+[dRes(3:2:end) 0]];
                         end
                         Resz=Resz+sparse(1,idx(logical(data.idx_perturb(:,j)),i),dRes,1,nz);
                    elseif j==2
                         if mod(data.nps,2)
                            dRes=dRes(1:2:end)+[dRes(2:2:end) 0];
                         else
                            dRes=dRes(1:2:end)+dRes(2:2:end);
                         end
                         Resz=Resz+sparse(1,idx(logical(data.idx_perturb(:,j)),i),dRes,1,nz);
                    else
                         Resz=Resz+sparse(1,idx(logical(data.idx_perturb(:,j)),i),dRes,1,nz);
                    end
                end
            else
                [ResNormCost_p,~]=costResidualMin_ModeMinRes( X,U,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)*T+t0+et0{i}*e,data);
                [ResNormCost_m,~]=costResidualMin_ModeMinRes( X,U,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)*T+t0-et0{i}*e,data);
                dRes=sum(ResNormCost_p-ResNormCost_m)/(2*e);
                Resz=Resz+sparse(1,idx(1,i),dRes,1,nz);
            end
        end
end

% Return the gradient
grad=Resz;

