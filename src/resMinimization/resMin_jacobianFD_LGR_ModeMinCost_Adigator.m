function jac=resMin_jacobianFD_LGR_ModeMinCost_Adigator(g,avrc,X_Np1,U_Np1,P,T_Np1,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,data,free_time,Res_vec)
%resMin_jacobianFD_LGR_ModeMinCost_Adigator - Jacobian computation for
%integrated residual minimization (alternating method: cost
%minimization) using Adigator
%
% Syntax:   jac=resMin_jacobianFD_LGR_ModeMinCost_Adigator(g,avrc,X_Np1,U_Np1,P,T_Np1,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,data,free_time,Res_vec)
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
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,~]=deal(dataNLP.sizes{1:17});
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


%% Compute contribution of R to the gradient
Resz=sparse(Res_vec.dY_location(:,1),Res_vec.dY_location(:,2),Res_vec.dY,Res_vec.dY_size(1),Res_vec.dY_size(2));


% Map derivatives to the jacobian
%---------------------------------
jac=[gz(dataNLP.gAllidx,:);rcz;bz;Resz];





