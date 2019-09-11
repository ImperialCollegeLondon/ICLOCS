function [ cr ] = addRateConstraint( x,u,p,t,data )
%addRateConstraint - Format the implementation for rate constraints.
%
% Syntax:  [ cr ] = addRateConstraint( x,u,p,t,data )
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    data- structured variable containing the values of additional data used inside
%          the function
% 
% Output:
%    cr - rate constraint function
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

if (strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))
    [nrcl,nrcu,nrce]=deal(data.sizes{14:16});
else
    [nrcl,nrcu,nrce]=deal(data.sizes{10:12});
end

nrc=nrcl+nrcu+nrce;

nx=size(x,2);
nu=size(u,2);

if nrc
    if strcmp(data.options.discretization,'hermite')
        ns=data.sizes{9};
        M=(size(x,1)-1)/ns;
        cu=zeros(M,nrcu/M);
        cl=zeros(M,nrcl/M);
        ce=zeros(M,nrce/M);
        DT=diff(t(1:2:end));
        k=1;l=1;m=1;
        for i=1:nx
            if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
                cl(:,k:k+2)=[data.RCmap.AxHS1*x(:,i)-DT.*data.data.xrl(i),data.RCmap.AxHS2*x(:,i)-DT.*data.data.xrl(i),data.RCmap.AxHS3*x(:,i)-DT.*data.data.xrl(i)];
                k=k+3;
            end
            if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
                ce(:,l:l+1)=[data.RCmap.AxHS3*x(:,i)-data.RCmap.DT.*data.data.xrl(i),data.RCmap.AxHS4*x(:,i)-0.5*DT.*data.data.xrl(i)];
                l=l+2;
            end
            if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
                cu(:,m:m+2)=[data.RCmap.AxHS1*x(:,i)-DT.*data.data.xru(i),data.RCmap.AxHS2*x(:,i)-DT.*data.data.xru(i),data.RCmap.AxHS3*x(:,i)-DT.*data.data.xru(i)];
                m=m+3;
            end
        end
        for i=1:nu
            if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
                cl(:,k:k+1)=[data.RCmap.AuHS1*u(:,i)-DT.*data.data.url(i),data.RCmap.AuHS2*u(:,i)-DT.*data.data.url(i)];
                k=k+2;
            end
            if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
                ce(:,l:l+1)=[data.RCmap.AuHS3*u(:,i)-DT.*data.data.url(i),data.RCmap.AuHS4*u(:,i)-0.5*DT.*data.data.url(i)];
                l=l+2;
            end
            if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
                cu(:,m:m+1)=[data.RCmap.AuHS1*u(:,i)-DT.*data.data.uru(i),data.RCmap.AuHS2*u(:,i)-DT.*data.data.uru(i)];
                m=m+2;
            end
        end

        cl=reshape(cl',size(cl,1)*size(cl,2),1)';
        cu=reshape(cu',size(cu,1)*size(cu,2),1)';
        ce=reshape(ce',size(ce,1)*size(ce,2),1)';


    elseif (strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))
        npd=data.sizes{10};
        M=sum(npd);
        t0=t(1);
        tf=t(end);
        AxLGR=data.map.D_structure;
        t_segment=(tf-t0)/2*data.tau_segment+(tf+t0)/2;
        DT_seg=zeros(1,M);
        DT_s=diff(t_segment);
        for i=1:length(npd)
            DT_seg(1+sum(npd(1:(i-1))):sum(npd(1:i)))=repmat(DT_s(i),1,npd(i));
        end
        DT_seg=DT_seg';
        
        cu=zeros(M,nrcu/M);
        cl=zeros(M,nrcl/M);
        ce=zeros(M,nrce/M);
        k=1;l=1;m=1;
        for i=1:nx
            if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
                cl(:,k)=AxLGR*x(:,i)-0.5*DT_seg.*data.data.xrl(i);
                k=k+1;
            end
            if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
                ce(:,l)=AxLGR*x(:,i)-0.5*DT_seg.*data.data.xrl(i);
                l=l+1;
            end
            if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
                cu(:,m)=AxLGR*x(:,i)-0.5*DT_seg.*data.data.xru(i);
                m=m+1;
            end
        end
        for i=1:nu
            if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
                cl(:,k)=AxLGR*[u(:,i);u(end,i)]-0.5*DT_seg.*data.data.url(i);
                k=k+1;
            end
            if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
                ce(:,l)=AxLGR*[u(:,i);u(end,i)]-0.5*DT_seg.*data.data.url(i);
                l=l+1;
            end
            if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
                cu(:,m)=AxLGR*[u(:,i);u(end,i)]-0.5*DT_seg.*data.data.uru(i);
                m=m+1;
            end
        end

        cl=cl(:)';
        cu=cu(:)';
        ce=ce(:)';

    else
        DT=diff(t);
        M=size(x,1)-1;
        cu=zeros(M,nrcu/M);
        cl=zeros(M,nrcl/M);
        ce=zeros(M,nrce/M);
        k=1;l=1;m=1;
        for i=1:nx
            if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
                cl(:,k)=data.RCmap.Ax*x(:,i)-DT.*data.data.xrl(i);
                k=k+1;
            end
            if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
                ce(:,l)=data.RCmap.Ax*x(:,i)-DT.*data.data.xrl(i);
                l=l+1;
            end
            if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
                cu(:,m)=data.RCmap.Ax*x(:,i)-DT.*data.data.xru(i);
                m=m+1;
            end
        end
        for i=1:nu
            if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
                cl(:,k)=data.RCmap.Au*u(:,i)-DT.*data.data.url(i);
                k=k+1;
            end
            if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
                ce(:,l)=data.RCmap.Au*u(:,i)-DT.*data.data.url(i);
                l=l+1;
            end
            if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
                cu(:,m)=data.RCmap.Au*u(:,i)-DT.*data.data.uru(i);
                m=m+1;
            end
        end
        cl=reshape(cl',size(cl,1)*size(cl,2),1)';
        cu=reshape(cu',size(cu,1)*size(cu,2),1)';
        ce=reshape(ce',size(ce,1)*size(ce,2),1)';

    end
    cr=[cl, cu, ce];
else
    cr=[];
end




end

