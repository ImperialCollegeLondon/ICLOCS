% -------------------------------------------------------------------------
%  2. Evaluate the Jacobian of the constraints
% -------------------------------------------------------------------------

function jac=jacobianFDAdigator_LGR(f,g,avrc,X_Np1,U,P,T,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data)
% jacobianFDAdigator_LGR - It evaluates the Jacobian of the constraints
% with Adigator
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
X=X_Np1(1:M,:);

% Compute fgz
%------------
fgdz=sparse(const_vec_Adigator.dY_location(:,1),const_vec_Adigator.dY_location(:,2),const_vec_Adigator.dY,const_vec_Adigator.dY_size(1),const_vec_Adigator.dY_size(2));

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
jac=[fgdz;rcz;bz];


