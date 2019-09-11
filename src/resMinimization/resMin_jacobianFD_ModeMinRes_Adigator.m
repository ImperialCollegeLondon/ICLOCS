function jac=resMin_jacobianFD_ModeMinRes_Adigator(L,E,g,avrc,X,Xr,U,Ur,P,T,b,x0,xf,u0,uf,p,t0,tf,data,Res_vec)
%resMin_jacobianFD_ModeMinRes_Adigator - Jacobian computation for
%integrated residual minimization (alternating method: residual
%minimization) using Adigator
%
% Syntax:   jac=resMin_jacobianFD_ModeMinRes_Adigator(L,E,g,avrc,X,Xr,U,Ur,P,T,b,x0,xf,u0,uf,p,t0,tf,data,Res_vec)
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
t=(tf-t0)*T+t0;

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
            rcp=avrc(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),data);
            rcm=avrc(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),data);
            rcz=rcz+sparse(1:nrc,idx(:,i),(rcp-rcm)/(2*e),nrc,nz);
        end
    end
    rcz=rcz+[dataNLP.map.Acl;dataNLP.map.Ae;dataNLP.map.Acu];
end

% Compute bz
%------------


if isfield(dataNLP.options,'resminRep') && dataNLP.options.resminRep.collmatch
    
    if nt
        bz=spalloc(2*n,nz,(2*m+2*n+nt+np)*2*n);
        idx=repmat(dataNLP.FD.index.Ey,2*n,1);nfd=size(idx,2);
        i_st=1;
        i_end=nfd;
    else
        bz=spalloc(2*n,nz,(2*m+2*n+np)*2*n);
        idx=repmat(dataNLP.FD.index.Ey,2*n,1);
        i_st=nt+1;
        i_end=nt+np+(m+n)*2;
        idx=idx-nt;
    end
    et0=e*dataNLP.FD.vector.Ey.et0;etf=e*dataNLP.FD.vector.Ey.etf;ep=e*dataNLP.FD.vector.Ey.ep;
    ex0=e*dataNLP.FD.vector.Ey.ex0;eu0=e*dataNLP.FD.vector.Ey.eu0;
    exf=e*dataNLP.FD.vector.Ey.exf;euf=e*dataNLP.FD.vector.Ey.euf;

    for i=i_st:i_end
        bp=collmatch_endpt(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(i),tf+etf(i),vdat);
        bm=collmatch_endpt(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0-et0(i),tf-etf(i),vdat);
        bz=bz+sparse(1:2*n,idx(:,i),(bp-bm)/(2*e),2*n,nz);
    end
    
    
else
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
end


% Compute contribution of L to the gradient
Lz=zeros(1,nz);

idx=dataNLP.FD.index.Ly;
nfd=size(idx,2);                               
et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;

if nt
    i_st=1;
    i_end=nfd;
else
    i_st=nt+1;
    i_end=nt+np+m+n;
    idx=idx-nt;
end

snm=ones(M,1);
for i=i_st:i_end
      dL=((tf+etf{i}*e-t0-et0{i}*e)*L(X+ex{i}*e,Xr,U+eu{i}*e,Ur,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)*T+t0+et0{i}*e,vdat)-...
      (tf-etf{i}*e-t0+et0{i}*e)*L(X-ex{i}*e,Xr,U-eu{i}*e,Ur,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)*T+t0-et0{i}*e,vdat))/(2*e);
      Lz=Lz+sparse(1,idx(:,i),dataNLP.map.W*dL,1,nz);
      JL{i}=dL;
end
if ~nfd
    JL=0;   
end   

% Compute contribution of E to the gradient
if nt
    Ez=spalloc(1,nz,nt+np+2*(n+m));
    idx=dataNLP.FD.index.Ey;
    nfd=size(idx,2);  
    i_st=1;
    i_end=nfd;
else
    Ez=spalloc(1,nz,np+2*(n+m));
    idx=dataNLP.FD.index.Ey;
    i_st=nt+1;
    i_end=nt+np+(m+n)*2;
    idx=idx-nt;
end

et0=e*dataNLP.FD.vector.Ey.et0;etf=e*dataNLP.FD.vector.Ey.etf;ep=e*dataNLP.FD.vector.Ey.ep;
ex0=e*dataNLP.FD.vector.Ey.ex0;eu0=e*dataNLP.FD.vector.Ey.eu0;
exf=e*dataNLP.FD.vector.Ey.exf;euf=e*dataNLP.FD.vector.Ey.euf;

for i=i_st:i_end
    Ez(1,idx(i))=(E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat)-...
              E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0-et0(:,i),tf-etf(:,i),vdat))/(2*e);
end



% Return the cost gradient
grad=Lz+Ez;


% Compute contribution of R to the gradient
Resz=sparse(Res_vec.dY_location(:,1),Res_vec.dY_location(:,2),Res_vec.dY,Res_vec.dY_size(1),Res_vec.dY_size(2));

% Map derivatives to the jacobian
%---------------------------------
jac=[gz(dataNLP.gAllidx,:);rcz;bz;grad;Resz];

