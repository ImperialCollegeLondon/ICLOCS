function [jac,Jf]=jacConstzLGR(df,dg,g,f,X,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data)

%jacConstzLGR - Return the Jacobian of the Constraints with respect to the 
%            variable z  when the analytic option has been selected
%
% Syntax:  [jac,Jf]=jacConstzLGR(df,dg,g,f,X,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data)
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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE ---------------


e=data.options.perturbation.J;                                % pertubation size
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{1:9});
nz=nt+np+(M+1)*n+M*m;
vdat=data.data;
% Compute fz
%------------

fz=spalloc(n*M,nz,n*M*n+m*M*n+np+nt);
Jaf=[];

if df.flag

  idx=data.FD.index.f;    
  % Compute and store the derivate of fz with respect to tf

  for i=1:n
    fz=fz+sparse(1:M*n,idx(:,i),reshape(repmat(data.t_segment_end,1,n).*df.dx{i},M*n,1),M*n,nz); 
  end
  for i=1:m
   fz=fz+sparse(1:M*n,idx(:,n+i),reshape(repmat(data.t_segment_end,1,n).*df.du{i},M*n,1),M*n,nz); 
  end
  Jaf=[df.dx, df.du];
  
  if np
      for i=1:np
        fz=fz+sparse(1:M*n,idx(:,n+m+i),reshape(repmat(data.t_segment_end,1,n).*df.dp{i},M*n,1),M*n,nz);
      end
  Jaf=[Jaf, df.dp];
  end
  
  if nt>=2 && ~isempty(df.dt)
     for i=1:nt
         if i==1
            fz_t=-f(X,U,P,(tf-t0)/2*T+(tf-t0)/2,vdat)*.5+repmat(data.t_segment_end,1,n).*(df.dt{i}).*(kron((1-T)/2,ones(1,n)));
         elseif i==nt
            fz_t=f(X,U,P,(tf-t0)/2*T+(tf-t0)/2,vdat)*.5+repmat(data.t_segment_end,1,n).*(df.dt{end}).*(kron((1+T)/2,ones(1,n)));
         end
         fz=fz+sparse(1:M*n,idx(:,n+m+np+i),reshape(fz_t.',M*n,1),M*n,nz);
     end
     Jaf=[Jaf, df.dt];
  end
  

  

  Jf=vertcat(Jaf);

else
%      idx=data.FD.index.f;nfd=size(idx,2);
%      etf=e*data.FD.vector.f.etf; et0=e*data.FD.vector.f.et0; ep=data.FD.vector.f.ep;et=e*data.FD.vector.f.et;
%      ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;
%  
%     for i=1:nfd
%        fp=f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
%        fm=f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
%        ft=((repmat(data.t_segment_end+data.t_segment_mat_m*(et{i}'*e),1,n)).*fp-(repmat(data.t_segment_end-data.t_segment_mat_m*(et{i}'*e),1,n)).*fm)/(2*e);
%        fz=fz+sparse(1:M*n,idx(:,i),reshape(ft,M*n,1),M*n,nz);
%        Jf{i}=(fp-fm)/(2*e);
%     end

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
            fp=(DTLP+etf(i)/2-et0(i)/2).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
            fm=(DTLP-etf(i)/2+et0(i)/2).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
            ft=sparse(1:M*n,idx(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
            fz=fz+ft;
        end

    end
end





% Compute gz
%------------

gz=spalloc(ng*M,nz,ng*M*(n+m+np+nt));
if ng
%     dg.flag=0;
   if dg.flag
    idx=data.FD.index.g;nfd=size(idx,2);
    dgz=[dg.dx dg.du];
    
      for i=1:(m+n)
        gz=gz+sparse(1:M*ng,idx(:,i),reshape(dgz{i},M*ng,1),M*ng,nz);
      end

      if np
          for i=1:np
            gz=gz+sparse(1:M*ng,idx(:,i+m+n),reshape(dg.dp{i},M*ng,1),M*ng,nz);  
          end
      end

      if nt>=2 && ~isempty(dg.dt)
         for i=1:nt
             if i==1
                 Tj=kron((1-T)/2,ones(1,ng)); 
             elseif i==nt
                 Tj=kron((T+1)/2,ones(1,ng)); 
             end
             
             gz=gz+sparse(1:M*ng,idx(:,i+m+n+np),reshape((dg.dt{i}.*Tj),M*ng,1),M*ng,nz);
         end
      end
  
    
   else    
%     idx=data.FD.index.g;nfd=size(idx,2);
%     et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
%     ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
%     for i=1:nfd
%       gp=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+(tf+etf(i)+t0+et0(i))/2,vdat);
%       gm=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+(tf-etf(i)+t0-et0(i))/2,vdat);
%       gz=gz+sparse(1:M*ng,idx(:,i),reshape((gp-gm)'/(2*e),M*ng,1),M*ng,nz);
%     end
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



% Compute bz
%------------

bz=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
if nb

    if db.flag 
        if data.options.adaptseg==1 
            bz(data.costStruct.B==1)=[db.dx0 db.dxf db.du0 db.duf db.dp db.dt];
        else
            bz(data.costStruct.B==1)=[db.dx0 db.dxf db.du0 db.duf db.dp db.dt0 db.dtf];
        end
    else     
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
end

% if ~isempty(find(bz))
%     bz
% end

% Map derivatives to the jacobian
%---------------------------------


% figure
% plot(fz)
% figure
% plot(gz)
% figure
% plot(bz)
% D_structure=kron(speye(nps),data.map.LGR.diff_matrix(:,1:end-1));
% for i=1:nps
%     D_structure((i-1)*npd+1:(i-1)*npd+npd,npd*i+1)=data.map.LGR.diff_matrix(:,end);
% end
% D_structure=[kron(speye(n*nps),data.map.LGR.diff_matrix(:,1:end-1)) zeros(M*n,M*m+np+nt)];
% repmat(data.t_segment_mat_m,n,1)
jac=[[kron(speye(n),data.map.D_structure) zeros(M*n,M*m+np+nt) ]-fz;gz(data.gAllidx,:);bz];






%------------- END CODE --------------