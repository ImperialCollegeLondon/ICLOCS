function jac=resMin_jacobianFD_ModeMinCost(g,avrc,X,U,P,T,b,x0,xf,u0,uf,p,t0,tf,data)
%resMin_jacobianFD_ModeMinCost - Jacobian computation for
%integrated residual minimization (alternating method: cost
%minimization) using finite difference
%
% Syntax:   jac=resMin_jacobianFD_ModeMinCost(g,avrc,X,U,P,T,b,x0,xf,u0,uf,p,t0,tf,data)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE ---------------


dataNLP=data.dataNLP;
e=dataNLP.options.perturbation.J;                                 % pertubation size
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,~]=deal(dataNLP.sizes{1:13});
nrc=nrcl+nrcu+nrce;
% 
if nt
    nz=nt+np+M*n+N*m;
else
    nz=np+M*n+N*m;
end
vdat=dataNLP.data;

% Compute gz
%------------
if nt
    gz=spalloc(ng*M,nz,ng*M*(n+m+np+nt));
else
    gz=spalloc(ng*M,nz,ng*M*(n+m+np));
end

if ng
    if nt
        idx=dataNLP.FD.index.g;
        nfd=size(idx,2);
        i_st=1;
        i_end=nfd;
    else
        idx=dataNLP.FD.index.g;
        i_st=nt+1;
        i_end=nt+np+m+n;
        idx=idx-nt;
    end
    et0=e*dataNLP.FD.vector.g.et0;etf=e*dataNLP.FD.vector.g.etf;ep=dataNLP.FD.vector.g.ep;
    ex=dataNLP.FD.vector.g.ex;eu=dataNLP.FD.vector.g.eu;

    for i=i_st:i_end

        gp=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),vdat);
        gm=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),vdat);
        gz=gz+sparse(1:M*ng,idx(:,i),reshape((gp-gm)'/(2*e),M*ng,1),M*ng,nz);

    end
end

% Compute rcz
%------------
if nt
    rcz=spalloc(nrc,nz,nrc*(n+m+np+nt));
else
    rcz=spalloc(nrc,nz,nrc*(n+m+np));
end
% 

if nrc
    if nt
        idx=dataNLP.FD.index.rc;nfd=size(idx,2);
        i_st=1;
        i_end=nfd;
    else
        idx=dataNLP.FD.index.rc;
        i_st=nt+1;
        i_end=nt+np+m+n;
        idx=idx-nt;
    end
    et0=e*dataNLP.FD.vector.rc.et0;etf=e*dataNLP.FD.vector.rc.etf;ep=dataNLP.FD.vector.rc.ep;
    ex=dataNLP.FD.vector.rc.ex;eu=dataNLP.FD.vector.rc.eu;
    
    for i=i_st:i_end
        if ~any(ex{i}(:)) && ~any(eu{i}(:))
            rcp=avrc(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),dataNLP);
            rcm=avrc(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),dataNLP);
            rcz=rcz+sparse(1:nrc,idx(:,i),(rcp-rcm)/(2*e),nrc,nz);
        end
    end
    rcz=rcz+[dataNLP.map.Acl;dataNLP.map.Ae;dataNLP.map.Acu];
end

% Compute bz
%------------
if nt
    bz=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
else
    bz=spalloc(nb,nz,(2*m+2*n+np)*nb);
end

if nb
    if nt
        idx=dataNLP.FD.index.b;nfd=size(idx,2);
        i_st=1;
        i_end=nfd;
    else
        idx=dataNLP.FD.index.b;
        i_st=nt+1;
        i_end=nt+np+(m+n)*2;
        idx=idx-nt;
    end
    et0=e*dataNLP.FD.vector.b.et0;etf=e*dataNLP.FD.vector.b.etf;ep=e*dataNLP.FD.vector.b.ep;
    ex0=e*dataNLP.FD.vector.b.ex0;eu0=e*dataNLP.FD.vector.b.eu0;
    exf=e*dataNLP.FD.vector.b.exf;euf=e*dataNLP.FD.vector.b.euf;

    for i=i_st:i_end
        bp=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(i),tf+etf(i),vdat);
        bm=b(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0-et0(i),tf-etf(i),vdat);
        bz=bz+sparse(1:nb,idx(:,i),(bp-bm)/(2*e),nb,nz);
    end
end


% Compute contribution of R to the gradient
idx=dataNLP.FD.index.Ly;
nfd=size(idx,2);                               
et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;
 
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
            [~,ResCost_p]=constResidualMin_ModeMinCost( X+ex{i}.*data.idx_perturb(:,j)*e,U+eu{i}.*data.idx_perturb(:,j)*e,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)*T+t0+et0{i}*e,data,1);
            [~,ResCost_m]=constResidualMin_ModeMinCost( X-ex{i}.*data.idx_perturb(:,j)*e,U-eu{i}.*data.idx_perturb(:,j)*e,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)*T+t0-et0{i}*e,data,1);
            dRes=(ResCost_p-ResCost_m)/(2*e);
            if j==1
                 if mod(data.nps,2)
                    dRes=[dRes(:,1) dRes(:,2:2:end)+dRes(:,3:2:end)];
                 else
                    dRes=[dRes(:,1) dRes(:,2:2:end)+[dRes(:,3:2:end) zeros(n,1)]];
                 end
                 idxl=idx(logical(data.idx_perturb(:,j)),i);
                 Resz=Resz+sparse(repmat(1:n,1,length(idxl)),repelem(idxl,n,1),dRes(:),n,nz);
            elseif j==2
                 if mod(data.nps,2)
                    dRes=dRes(:,1:2:end)+[dRes(:,2:2:end) zeros(n,1)];
                 else
                    dRes=dRes(:,1:2:end)+dRes(:,2:2:end);
                 end
                 idxl=idx(logical(data.idx_perturb(:,j)),i);
                 Resz=Resz+sparse(repmat(1:n,1,length(idxl)),repelem(idxl,n,1),dRes(:),n,nz);
            else
                 idxl=idx(logical(data.idx_perturb(:,j)),i);
                 Resz=Resz+sparse(repmat(1:n,1,length(idxl)),repelem(idxl,n,1),dRes(:),n,nz);
            end
        end
    else
        [~,ResCost_p]=constResidualMin_ModeMinCost( X,U,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)*T+t0+et0{i}*e,data,1);
        [~,ResCost_m]=constResidualMin_ModeMinCost( X,U,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)*T+t0-et0{i}*e,data,1);
        dRes=sum((ResCost_p-ResCost_m)/(2*e),2);
        idxl=idx(1,i);
        Resz=Resz+sparse(repmat(1:n,1,length(idxl)),repelem(idxl,n,1),dRes(:),n,nz);
    end
end


% Map derivatives to the jacobian
%---------------------------------
jac=[gz(dataNLP.gAllidx,:);rcz;bz;Resz];

