function [jac,Jf]=jacConstz(df,dg,g,f,avrc,X,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data)
%JACCONSTZ - Return the Jacobian of the Constraints with respect to the 
%            variable z  when the analytic option has been selected
%
% Syntax:  [jac,Jf]=jacConstz(df,dg,g,f,avrc,X,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data)
%
% Inputs:
%    see directCollocation.m
%
% Outputs:
%    jac, Jf - jacobian information
%
% Other m-files required: none
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

%------------- BEGIN CODE ---------------


e=data.options.perturbation.J;                                % pertubation size
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,~]=deal(data.sizes{1:13});
nrc=nrcl+nrcu+nrce;
nz=nt+np+M*n+N*m;
vdat=data.data;
% t0=data.t0;
% k0=data.k0;
% Compute fz
%------------

fz=spalloc(n*M,nz,n*M*(n+m)+np+nt);
Jaf=[];


if df.flag

  idx=data.FD.index.f;    
  % Compute and store the derivate of fz with respect to tf

   if nt==1 && ~isempty(df.dt) && ~isempty(cell2mat(df.dt))
     fz_tf=f(X,U,P,(tf-t0)*T+t0,vdat)+(tf-t0)*(df.dt{1}).*(kron(T,ones(1,n)));
     fz=sparse(1:M*n,idx(:,1),reshape(fz_tf.',M*n,1),M*n,nz);
     Jaf=[df.dt];
   elseif nt==2 && ~isempty(df.dt) && ~isempty(cell2mat(df.dt))
     fz_t0=-(f(X,U,P,(tf-t0)*T+t0,vdat)+(tf-t0)*(df.dt{1}).*(kron(T,ones(1,n))));
     fz_tf=f(X,U,P,(tf-t0)*T+t0,vdat)+(tf-t0)*(df.dt{1}).*(kron(T,ones(1,n)));
     fz=fz+sparse(1:M*n,idx(:,1),reshape(fz_t0.',M*n,1),M*n,nz)+sparse(1:M*n,idx(:,2),reshape(fz_tf.',M*n,1),M*n,nz);
     Jaf=[df.dt, df.dt];
   else
         idx=data.FD.index.f;
         et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
         ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

        for i=1:nt
           fp=f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),vdat);
           fm=f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),vdat);
           ft=((tf+etf(i)-t0-et0(i))*fp-(tf-etf(i)-t0+et0(i))*fm)/(2*e);
           fz=fz+sparse(1:M*n,idx(:,i),reshape(ft',M*n,1),M*n,nz);
        end
        Jaf=[cell(1,1)];
   end
  if np
  for i=1:np
    fz=fz+sparse(1:M*n,idx(:,nt+i),reshape((tf-t0)*df.dp{i},M*n,1),M*n,nz);
  end
   Jaf=[Jaf, df.dp];
  end
  for i=1:n
     
  fz=fz+sparse(1:M*n,idx(:,nt+np+i),reshape((tf-t0)*df.dx{i}',M*n,1),M*n,nz); 
  
  end
  for i=1:m
   fz=fz+sparse(1:M*n,idx(:,nt+np+n+i),reshape((tf-t0)*df.du{i}',M*n,1),M*n,nz); 
  end
  Jaf=[Jaf, df.dx, df.du];
  Jf=vertcat(Jaf);



else
     idx=data.FD.index.f;nfd=size(idx,2);
     et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
     ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    for i=1:nfd
       fp=f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),vdat);
       fm=f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),vdat);
       ft=((tf+etf(i)-t0-et0(i))*fp-(tf-etf(i)-t0+et0(i))*fm)/(2*e);
       fz=fz+sparse(1:M*n,idx(:,i),reshape(ft',M*n,1),M*n,nz);
       Jf{i}=(fp-fm)/(2*e);
    end

end



% Compute gz
%------------

gz=spalloc(ng*M,nz,ng*M*(n+m+np+nt));
if ng
   if dg.flag
    idx=data.FD.index.g;nfd=size(idx,2);
    dgz=[dg.dx dg.du];
    for i=1:nfd
      if nt&&(i<=nt)
        Tj=kron(T,ones(1,ng)); 
        gz=gz+sparse(1:M*ng,idx(:,i),reshape((dg.dt{i}.*Tj)',M*ng,1),M*ng,nz);
      elseif np&&(nt<i)&&(i<nt+np+1)
        gz=gz+sparse(1:M*ng,idx(:,i),reshape(dg.dp{i-nt}',M*ng,1),M*ng,nz);  
      else   
       gz=gz+sparse(1:M*ng,idx(:,i),reshape(dgz{i-nt-np}',M*ng,1),M*ng,nz);
      end
    end
    
   else    
    idx=data.FD.index.g;nfd=size(idx,2);
    et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
    ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
    for i=1:nfd
      gp=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),vdat);
      gm=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),vdat);
      gz=gz+sparse(1:M*ng,idx(:,i),reshape((gp-gm)'/(2*e),M*ng,1),M*ng,nz);
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
            rcp=avrc(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),data);
            rcm=avrc(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),data);
            rcz=rcz+sparse(1:nrc,idx(:,i),(rcp-rcm)/(2*e),nrc,nz);
        end
    end
    rcz=rcz+[data.map.Acl;data.map.Ae;data.map.Acu];
end

% Compute bz
%------------

bz=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
if nb

if db.flag 
    if nt==1
        bz(data.costStruct.B==1)=[db.dtf db.dp db.dx0 db.du0 db.duf db.dxf];
    elseif nt==2
        bz(data.costStruct.B==1)=[db.dt0 db.dtf db.dp db.dx0 db.du0 db.duf db.dxf];
    else
        bz(data.costStruct.B==1)=[db.dp db.dx0 db.du0 db.duf db.dxf];
    end
 else     
 idx=data.FD.index.b;nfd=size(idx,2);
 et0=e*data.FD.vector.b.et0;etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;
 ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
 exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;
 for i=1:nfd
    bp=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0,tf+etf,vdat);
    bm=b(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0-et0,tf-etf,vdat);
    bz=bz+sparse(1:nb,idx(:,i),(bp-bm)/(2*e),nb,nz);
 end
end
end

% Map derivatives to the jacobian
%---------------------------------
jac=[[zeros(n,nt) zeros(n,np) eye(n), zeros(n,(M-1)*n+N*m)]*data.cx0;data.map.A*data.map.Vx+data.map.B*fz;gz;rcz;bz];






%------------- END CODE --------------