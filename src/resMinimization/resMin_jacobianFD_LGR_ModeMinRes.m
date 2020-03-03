function jac=resMin_jacobianFD_LGR_ModeMinRes(L,E,g,avrc,X_Np1,Xr,U_Np1,Ur,P,T_Np1,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_mat_m,data,free_time)
%resMin_jacobianFD_LGR_ModeMinRes - Jacobian computation for
%integrated residual minimization (alternating method: residual
%minimization) using finite difference
%
% Syntax:   jac=resMin_jacobianFD_LGR_ModeMinRes(L,E,g,avrc,X_Np1,Xr,U_Np1,Ur,P,T_Np1,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_mat_m,data,free_time)
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
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,~,ng_eq,ng_neq]=deal(dataNLP.sizes{1:19});
nrc=nrcl+nrcu+nrce;
nx=(M+1)*n;                               % Number of unknown states
nu=M*m;                               % Number of unknown controls
% Length of the optimization variable 
if free_time
    nz=nt+np+nx+nu;
else
    nz=nx+nu+np;
end

vdat=dataNLP.data;
DTLP=repmat(t_segment_end,1,n);

X=X_Np1(1:M,:);
U=U_Np1(1:M,:);
T=T_Np1(1:M,:);

% Compute gz
%------------
gz=spalloc(ng*M,nz,ng*M*(n+m+np+nt));
if ng
    idx=dataNLP.FD.index.g;nfd=size(idx,2);
    et0=e*dataNLP.FD.vector.g.et0;etf=e*dataNLP.FD.vector.g.etf;ep=dataNLP.FD.vector.g.ep;
    ex=dataNLP.FD.vector.g.ex;eu=dataNLP.FD.vector.g.eu;
    
    if free_time
        i_st=1;
        i_end=nfd;
    else
        i_st=1;
        i_end=m+n+np;
    end
       
    for i=i_st:i_end
        gp=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
        gm=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
        gz=gz+sparse(1:M*ng,idx(:,i),reshape((gp-gm)/(2*e),M*ng,1),M*ng,nz);
    end
end


% Compute rcz
%------------
rcz=spalloc(nrc,nz,nrc*(n+m+np+nt));
if nrc

    idx=dataNLP.FD.index.rc;nfd=size(idx,2);
    et0=e*dataNLP.FD.vector.rc.et0;etf=e*dataNLP.FD.vector.rc.etf;ep=dataNLP.FD.vector.rc.ep;
    ex=dataNLP.FD.vector.rc.ex;eu=dataNLP.FD.vector.rc.eu;

   
    if free_time
        i_st=1;
        i_end=nfd;
    else
        i_st=1;
        i_end=m+n+np;
    end
    for i=i_st:i_end
        if ~any(ex{i}(:)) && ~any(eu{i}(:))
            rcp=avrc(X_Np1+ex{i}*e,U+eu{i}*e,P+ep{i}*e,[(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2;tf],dataNLP);
            rcm=avrc(X_Np1-ex{i}*e,U-eu{i}*e,P-ep{i}*e,[(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2;tf],dataNLP);
            rcz=rcz+sparse(1:nrc,idx(:,i),(rcp-rcm)/(2*e),nrc,nz);
        end
    end
    rcz=rcz+[dataNLP.map.Acl(:,1:end-nt);dataNLP.map.Ae(:,1:end-nt);dataNLP.map.Acu(:,1:end-nt)];
end


% Compute bz
%------------
if isfield(dataNLP.options,'resminRep') && dataNLP.options.resminRep.collmatch
    
    bz=spalloc(2*n,nz,(2*m+2*n+nt+np)*2*n);
    idx=repmat(dataNLP.FD.index.Ey,2*n,1);nfd=size(idx,2);
    et0=e*dataNLP.FD.vector.Ey.et0;etf=e*dataNLP.FD.vector.Ey.etf;ep=e*dataNLP.FD.vector.Ey.ep;
    ex0=e*dataNLP.FD.vector.Ey.ex0;eu0=e*dataNLP.FD.vector.Ey.eu0;
    exf=e*dataNLP.FD.vector.Ey.exf;euf=e*dataNLP.FD.vector.Ey.euf;

    if free_time
        i_st=1;
        i_end=nfd;
    else
        i_st=1;
        i_end=(m+n)*2;
    end

    if dataNLP.options.adaptseg==1 
        et=e.*dataNLP.FD.vector.Ey.et;
        for i=i_st:i_end
            t_segment_p=dataNLP.t_segment+et(i,:)';
            t_segment_m=dataNLP.t_segment-et(i,:)';
            bp=collmatch_endpt(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat,dataNLP.options,t_segment_p);
            bm=collmatch_endpt(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat,dataNLP.options,t_segment_m);
            bz=bz+sparse(1:2*n,idx(:,i),(bp-bm)/(2*e),2*n,nz);
        end
    else
        for i=i_st:i_end
        bp=collmatch_endpt(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat,dataNLP.options,[]);
        bm=collmatch_endpt(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat,dataNLP.options,[]);
        bz=bz+sparse(1:2*n,idx(:,i),(bp-bm)/(2*e),2*n,nz);
        end
    end
    
    
else

    bz=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
    if nb

    idx=dataNLP.FD.index.b;nfd=size(idx,2);
    et0=e*dataNLP.FD.vector.b.et0;etf=e*dataNLP.FD.vector.b.etf;ep=e*dataNLP.FD.vector.b.ep;et=e.*dataNLP.FD.vector.b.et;
    ex0=e*dataNLP.FD.vector.b.ex0;eu0=e*dataNLP.FD.vector.b.eu0;
    exf=e*dataNLP.FD.vector.b.exf;euf=e*dataNLP.FD.vector.b.euf;

    if free_time
        i_st=1;
        i_end=nfd;
    else
        i_st=1;
        i_end=(m+n)*2;
    end

    if dataNLP.options.adaptseg==1 
        for i=i_st:i_end
            t_segment_p=dataNLP.t_segment+et(i,:)';
            t_segment_m=dataNLP.t_segment-et(i,:)';
            bp=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat,dataNLP.options,t_segment_p);
            bm=b(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat,dataNLP.options,t_segment_m);
            bz=bz+sparse(1:nb,idx(:,i),(bp-bm)/(2*e),nb,nz);
        end
    else
        for i=i_st:i_end
        bp=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat,dataNLP.options,[]);
        bm=b(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat,dataNLP.options,[]);
        bz=bz+sparse(1:nb,idx(:,i),(bp-bm)/(2*e),nb,nz);
        end
    end

    end
end


% Compute contribution of L to the gradient

Lz=zeros(1,nz);

if dataNLP.FD.index.dL.flag==1          %  dataNLP.FD.index.dE.flag==1 when the analytic expression  of the gradient for the stage cost L is supplied  
    [dL,dE]=gradCost(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,tf,dataNLP);
    
    idx=dataNLP.FD.index.Ly;
    for  i=1:n
       dLx=dataNLP.map.w.*(dataNLP.t_segment_end.*dL.dx(:,i));
       Lz=Lz+sparse(1,idx(:,i),dLx,1,nz);
       JaL{i}=dL.dx(:,i);
    end  
    

    for i=1:m
       dLu=dataNLP.map.w.*(dataNLP.t_segment_end.*dL.du(:,i));  
       Lz=Lz+sparse(1,idx(:,n+i),dLu,1,nz);
       JaL{i+n}=dL.du(:,i);
    end
    

    if np
        for i=1:np
            dLp=dataNLP.map.w*(dataNLP.t_segment_end.*dL.dp(:,i));
            Lz=Lz+sparse(1,idx(:,n+m+i),dLp,1,nz);
            JaL{i+n+m}=dL.dp(:,i);
        end
    end
    

    if nt==2
        for i=1:nt
            if i==1
                dLt=-0.5*dataNLP.map.w'*L(X,Xr,U,Ur,P,t,vdat)+dataNLP.map.w'*(dataNLP.t_segment_end.*(dL.dt.*(1-T)/2)); %t0
            elseif i==2
                dLt=0.5*dataNLP.map.w'*L(X,Xr,U,Ur,P,t,vdat)+dataNLP.map.w'*(dataNLP.t_segment_end.*(dL.dt.*(1+T)/2)); %tf
            end
            Lz=Lz+sparse(1,idx(:,n+m+np+i),dLt,1,nz);
            JaL{i+n+m+np}=dL.dt(:);
        end
    end
   
    
    JL=vertcat(JaL);

else   % Numerical evalution of the gradient of the stage cost
     idx=dataNLP.FD.index.Ly;
     nfd=size(idx,2);   
     et=dataNLP.FD.vector.Ly.et;et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;
     ep=dataNLP.FD.vector.Ly.ep;
     ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;

    if free_time
        i_st=1;
        i_end=nfd;
    else
        i_st=1;
        i_end=m+n+np;
    end

     %This vector is used to adjust the size of the vector of the derivative
     %when the stage cost is identically zero 

     if dataNLP.options.adaptseg==1 
         snm=ones(M,1);
         for i=i_st:i_end
          dL=(((t_segment_end+t_segment_mat_m*(et{i}'*e)).*L(X+ex{i}*e,Xr,U+eu{i}*e,Ur,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*T+(tf+etf{i}*e+t0+et0{i}*e)/2,vdat)-...
          (t_segment_end-t_segment_mat_m*(et{i}'*e)).*L(X-ex{i}*e,Xr,U-eu{i}*e,Ur,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*T+(tf-etf{i}*e+t0-et0{i}*e)/2,vdat)).*snm)/(2*e);
          Lz=Lz+sparse(1,idx(:,i),diag(dataNLP.map.w)*dL,1,nz);
          JL{i}=dL;
         end
     else
         snm=ones(M,1);
         for i=i_st:i_end
          dL=(((t_segment_end+etf{i}*e/2-et0{i}*e/2).*L(X+ex{i}*e,Xr,U+eu{i}*e,Ur,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*T+(tf+etf{i}*e+t0+et0{i}*e)/2,vdat)-...
          (t_segment_end-etf{i}*e/2+et0{i}*e/2).*L(X-ex{i}*e,Xr,U-eu{i}*e,Ur,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*T+(tf-etf{i}*e+t0-et0{i}*e)/2,vdat)).*snm)/(2*e);
          Lz=Lz+sparse(1,idx(:,i),dataNLP.map.w.*dL,1,nz);
          JL{i}=dL;
         end
     end

    if ~nfd
     JL=0;  
    end   
end



% Compute contribution of E to the gradient
if free_time
    Ez=spalloc(1,nz,nt+np+2*(n+m));
    idx=dataNLP.FD.index.Ey;
    nfd=size(idx,2); 
    i_st=1;
    i_end=nfd;
else
    Ez=spalloc(1,nz,2*(n+m)+np);
    idx=dataNLP.FD.index.Ey;
    nfd=size(idx,2); 
    i_st=1;
    i_end=(m+n)*2+np;
end

if dataNLP.FD.index.dE.flag==1
    [dL,dE]=gradCost(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,tf,dataNLP);    
    Ez(dataNLP.costStruct.E==1)=[dE.dx0 dE.dxf dE.du0 dE.duf dE.dp dE.dt0 dE.dtf ];

else    
                                  
    et0=e*dataNLP.FD.vector.Ey.et0;etf=e*dataNLP.FD.vector.Ey.etf;ep=e*dataNLP.FD.vector.Ey.ep;
    ex0=e*dataNLP.FD.vector.Ey.ex0;eu0=e*dataNLP.FD.vector.Ey.eu0;
    exf=e*dataNLP.FD.vector.Ey.exf;euf=e*dataNLP.FD.vector.Ey.euf;

    for i=i_st:i_end
    Ez(1,idx(i))=(E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat)-...
                  E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat))/(2*e);
    end

   
end

% Return the gradient
grad=Lz+Ez;

%% Compute contribution of R to the gradient

idx=dataNLP.FD.index.Ly;
nfd=size(idx,2);                               
etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;

nx=(M+1)*n;                               % Number of unknown states
nu=M*m;                               % Number of unknown controls
et0=dataNLP.FD.vector.Ly.et0;


if data.free_time
    i_st=1;
    i_end=nfd;
    nz=nx+nu+nt+np;                       % Length of the optimization variable 
else
    i_st=1;
    i_end=m+n+np;
    nz=nx+nu+np;                       % Length of the optimization variable 
end

Resz=sparse(n+ng_eq,nz);
tau_seg_idx=[dataNLP.tau_seg_idx;dataNLP.tau_seg_idx(end)];
idx_state=[idx(:,1:n);idx(end,1:n)+1];
for i=i_st:i_end
    if (data.free_time && i<=(m+n)) || (~data.free_time)
        for j=1:size(data.idx_perturb,2)
            if i<=n
                [~,ResCost_p]=costResidualMin_ModeMinRes( X_Np1+[ex{i};ex{i}(end,:)].*data.idx_perturb(:,j)*e,U_Np1+[eu{i};eu{i}(end,:)].*data.idx_perturb(:,j)*e,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*[T;1]+(tf+etf{i}*e+t0+et0{i}*e)/2,data);
                [~,ResCost_m]=costResidualMin_ModeMinRes( X_Np1-[ex{i};ex{i}(end,:)].*data.idx_perturb(:,j)*e,U_Np1-[eu{i};eu{i}(end,:)].*data.idx_perturb(:,j)*e,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*[T;1]+(tf-etf{i}*e+t0-et0{i}*e)/2,data);
                dRes=(ResCost_p-ResCost_m)/(2*e);
            else
                [~,ResCost_p]=costResidualMin_ModeMinRes( X_Np1+[ex{i};ex{i}(end,:)].*data.idx_perturb(:,j)*e,U_Np1+[eu{i};eu{i}(end,:)].*[data.idx_perturb(1:end-1,j);data.idx_perturb(end-1,j)]*e,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*[T;1]+(tf+etf{i}*e+t0+et0{i}*e)/2,data);
                [~,ResCost_m]=costResidualMin_ModeMinRes( X_Np1-[ex{i};ex{i}(end,:)].*data.idx_perturb(:,j)*e,U_Np1-[eu{i};eu{i}(end,:)].*[data.idx_perturb(1:end-1,j);data.idx_perturb(end-1,j)]*e,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*[T;1]+(tf-etf{i}*e+t0-et0{i}*e)/2,data);
                dRes=(ResCost_p-ResCost_m)/(2*e);
            end
            if j==1
                 if mod(data.nps,2)
                    dRes=[dRes(:,1) dRes(:,2:2:end)+dRes(:,3:2:end)];
                 else
                    dRes=[dRes(:,1) dRes(:,2:2:end)+[dRes(:,3:2:end) zeros(n+ng_eq,1)]];
                 end
                 if i<=n
                    idxl=idx_state(logical(data.idx_perturb(:,j)),i);
                    Resz=Resz+sparse(repmat(1:n+ng_eq,1,length(idxl)),repelem(idxl,n+ng_eq,1),dRes(:),n+ng_eq,nz);
                 else
                    idxl=idx(logical(data.idx_perturb(1:end-1,j)),i);
                    dResl=dRes(:,1:length(idxl));
                    Resz=Resz+sparse(repmat(1:n+ng_eq,1,length(idxl)),repelem(idxl,n+ng_eq,1),dResl(:),n+ng_eq,nz);
                 end
            elseif j==2
                 if mod(data.nps,2)
                    dRes=dRes(:,1:2:end)+[dRes(:,2:2:end) zeros(n+ng_eqn,1)];
                 else
                    dRes=dRes(:,1:2:end)+dRes(:,2:2:end);
                 end
                 if i<=n
                    idxl=idx_state(logical(data.idx_perturb(:,j)),i);
                    Resz=Resz+sparse(repmat(1:n+ng_eq,1,length(idxl)),repelem(idxl,n+ng_eq,1),dRes(:),n+ng_eq,nz);
                 else
                    idxl=idx(logical(data.idx_perturb(1:end-1,j)),i);
                    dResl=dRes(:,1:length(idxl));
                    Resz=Resz+sparse(repmat(1:n+ng_eq,1,length(idxl)),repelem(idxl,n+ng_eq,1),dResl(:),n+ng_eq,nz);
                 end
            else
                 dRes=dRes(:,tau_seg_idx(logical(data.idx_perturb(:,j))));
                 if i<=n
                    idxl=idx_state(logical(data.idx_perturb(:,j)),i);
                    Resz=Resz+sparse(repmat(1:n+ng_eq,1,length(idxl)),repelem(idxl,n+ng_eq,1),dRes(:),n+ng_eq,nz);
                 else
                    idxl=idx(logical(data.idx_perturb(1:end-1,j)),i);
                    dResl=dRes(:,1:length(idxl));
                    Resz=Resz+sparse(repmat(1:n+ng_eq,1,length(idxl)),repelem(idxl,n+ng_eq,1),dResl(:),n+ng_eq,nz);
                 end
            end
        end
    else
        [~,ResCost_p]=costResidualMin_ModeMinRes( X_Np1,U_Np1,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*[T;1]+(tf+etf{i}*e+t0+et0{i}*e)/2,data);
        [~,ResCost_m]=costResidualMin_ModeMinRes( X_Np1,U_Np1,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*[T;1]+(tf-etf{i}*e+t0-et0{i}*e)/2,data);
        dRes=sum((ResCost_p-ResCost_m)/(2*e),2);
        idxl=idx(1,i);
        Resz=Resz+sparse(repmat(1:n+ng_eq,1,length(idxl)),repelem(idxl,n+ng_eq,1),dRes(:),n+ng_eq,nz);
    end
end



% Map derivatives to the jacobian
%---------------------------------
jac=[gz(dataNLP.gAllidx,:);rcz;bz;grad;Resz];





