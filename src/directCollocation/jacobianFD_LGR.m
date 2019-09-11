function jac=jacobianFD_LGR(f,g,avrc,X_Np1,U,P,T,b,x0,xf,u0,uf,p,t0,tf,data)

% jacobianFD_LGR - It evaluates numerically the Jacobian of the constraints
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk



e=data.options.perturbation.J;                                 % pertubation size
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,~]=deal(data.sizes{1:17});
nrc=nrcl+nrcu+nrce;
nz=nt+np+(M+1)*n+M*m;
vdat=data.data;
DTLP=repmat(data.t_segment_end,1,n);
DT_ratio_diff=repmat(data.tau_segment_ratio_diff,1,n);
X=X_Np1(1:M,:);
fg=vdat.functionfg;

% Compute fz and gz
%------------

fz=spalloc(n*M,nz,n*M*(n+m)+np+nt);
gz=spalloc(ng*M,nz,ng*M*(n+m+np+nt));

if ng && size(data.FD.index.f,2)==size(data.FD.index.g,2)
    idxf=data.FD.index.f;idxg=data.FD.index.g;nfd=size(idxf,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;et=data.FD.vector.f.et;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    
    for i=1:nfd
        [dyn_p,gp]=fg(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
        [dyn_m,gm]=fg(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
        if data.options.adaptseg==1 
            fp=(repmat(data.t_segment_end+data.t_segment_mat_m*(et{i}'*e),1,n)).*dyn_p;
            fm=(repmat(data.t_segment_end-data.t_segment_mat_m*(et{i}'*e),1,n)).*dyn_m;
            fz=fz+sparse(1:M*n,idxf(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
        else
            fp=(DTLP+etf(i)/2*DT_ratio_diff-et0(i)/2*DT_ratio_diff).*dyn_p;
            fm=(DTLP-etf(i)/2*DT_ratio_diff+et0(i)/2*DT_ratio_diff).*dyn_m;
            ft=sparse(1:M*n,idxf(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
            fz=fz+ft;
        end
        gz=gz+sparse(1:M*ng,idxg(:,i),reshape((gp-gm)/(2*e),M*ng,1),M*ng,nz);
    end
    
else
    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;et=data.FD.vector.f.et;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    if data.options.adaptseg==1 
        for i=1:nfd
            fp=(repmat(data.t_segment_end+data.t_segment_mat_m*(et{i}'*e),1,n)).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
            fm=(repmat(data.t_segment_end-data.t_segment_mat_m*(et{i}'*e),1,n)).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
            fz=fz+sparse(1:M*n,idx(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
        end
    else
        for i=1:nfd
            fp=(DTLP+etf(i)/2*DT_ratio_diff-et0(i)/2*DT_ratio_diff).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
            fm=(DTLP-etf(i)/2*DT_ratio_diff+et0(i)/2*DT_ratio_diff).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
            ft=sparse(1:M*n,idx(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
            fz=fz+ft;
        end
    end

    if ng
        idx=data.FD.index.g;nfd=size(idx,2);
        et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
        ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
        for i=1:nfd
            gp=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
            gm=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
            gz=gz+sparse(1:M*ng,idx(:,i),reshape((gp-gm)/(2*e),M*ng,1),M*ng,nz);
        end
    end
end

% Compute rcz
%------------

rcz=spalloc(nrc,nz,nrc*(n+m+np+nt));
if nrc

    idx=data.FD.index.rc;nfd=size(idx,2);
    et0=e*data.FD.vector.rc.et0;etf=e*data.FD.vector.rc.etf;ep=data.FD.vector.rc.ep;
    ex=data.FD.vector.rc.ex;eu=data.FD.vector.rc.eu;

   
for i=1:nfd
    if ~any(ex{i}(:)) && ~any(eu{i}(:))
        rcp=avrc(X_Np1+ex{i}*e,U+eu{i}*e,P+ep{i}*e,[(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2;tf],data);
        rcm=avrc(X_Np1-ex{i}*e,U-eu{i}*e,P-ep{i}*e,[(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2;tf],data);
        rcz=rcz+sparse(1:nrc,idx(:,i),(rcp-rcm)/(2*e),nrc,nz);
    end
end
    rcz=rcz+[data.map.Acl;data.map.Ae;data.map.Acu];
end


% Compute bz
%------------
bz=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
if nb

idx=data.FD.index.b;nfd=size(idx,2);
et0=e*data.FD.vector.b.et0;etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;et=e.*data.FD.vector.b.et;
ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;

if data.options.adaptseg==1 
    for i=1:nfd
        t_segment_p=data.t_segment+et(i,:)';
        t_segment_m=data.t_segment-et(i,:)';
        bp=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat,data.options,t_segment_p);
        bm=b(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat,data.options,t_segment_m);
        bz=bz+sparse(1:nb,idx(:,i),(bp-bm)/(2*e),nb,nz);
    end
else
    for i=1:nfd
    bp=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat,data.options,[]);
    bm=b(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat,data.options,[]);
    bz=bz+sparse(1:nb,idx(:,i),(bp-bm)/(2*e),nb,nz);
    end
end

end

% Map derivatives to the jacobian
%---------------------------------
jac=[[kron(speye(n),data.map.D_structure) zeros(M*n,M*m+np+nt)]+fz;gz(data.gAllidx,:);rcz;bz];





