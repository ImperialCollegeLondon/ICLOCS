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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

nx=size(x,2);
nu=size(u,2);
M=size(x,1);


cu=[];
cl=[];
ce=[];
if strcmp(data.options.transcription,'hermite')

    DT=diff(t(1:2:end));
    N=size(u,1);
    
    AxHS1=diag(-3*ones(M,1))+diag(4*ones(M-1,1),1)+diag(-ones(M-2,1),2);
    AxHS1([2:2:end,end],:)=[];
    AxHS2=diag(ones(M,1))+diag(-4*ones(M-1,1),1)+diag(3*ones(M-2,1),2);
    AxHS2([2:2:end,end],:)=[];
    AxHS3=diag(-ones(M,1))+diag(ones(M-2,1),2);
    AxHS3([2:2:end,end],:)=[];
    AxHS4=diag(-ones(M,1))+diag(ones(M-1,1),1);
    AxHS4([2:2:end,end],:)=[];
    
    AuHS1=diag(-3*ones(N,1))+diag(4*ones(N-1,1),1)+diag(-ones(N-2,1),2);
    AuHS1([2:2:end,end],:)=[];
    AuHS2=diag(ones(N,1))+diag(-4*ones(N-1,1),1)+diag(3*ones(N-2,1),2);
    AuHS2([2:2:end,end],:)=[];
    AuHS3=diag(-ones(N,1))+diag(ones(N-2,1),2);
    AuHS3([2:2:end,end],:)=[];
    AuHS4=diag(-ones(N,1))+diag(ones(N-1,1),1);
    AuHS4([2:2:end,end],:)=[];
    

    for i=1:nx
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            cl=[cl,[AxHS1*x(:,i)-DT.*data.data.xrl(i),AxHS2*x(:,i)-DT.*data.data.xrl(i),AxHS3*x(:,i)-DT.*data.data.xrl(i)]];
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            ce=[ce,[AxHS3*x(:,i)-DT.*data.data.xrl(i),AxHS4*x(:,i)-0.5*DT.*data.data.xrl(i)]];
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            cu=[cu,[AxHS1*x(:,i)-DT.*data.data.xru(i),AxHS2*x(:,i)-DT.*data.data.xru(i),AxHS3*x(:,i)-DT.*data.data.xru(i)]];
        end
    end
    for i=1:nu
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            cl=[cl,[AuHS1*u(:,i)-DT.*data.data.url(i),AuHS2*u(:,i)-DT.*data.data.url(i)]];
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            ce=[ce,[AuHS3*u(:,i)-DT.*data.data.url(i),AuHS4*u(:,i)-0.5*DT.*data.data.url(i)]];
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            cu=[cu,[AuHS1*u(:,i)-DT.*data.data.uru(i),AuHS2*u(:,i)-DT.*data.data.uru(i)]];
        end
    end

    cl=reshape(cl',size(cl,1)*size(cl,2),1)';
    cu=reshape(cu',size(cu,1)*size(cu,2),1)';
    ce=reshape(ce',size(ce,1)*size(ce,2),1)';


elseif (strcmp(data.options.transcription,'globalLGR')) || (strcmp(data.options.transcription,'hpLGR'))
    npd=data.sizes{10};
    t0=t(1);
    tf=t(end);
    AxLGR=data.map.D_structure;
    t_segment=(tf-t0)/2*data.tau_segment+(tf+t0)/2;
    DT_seg=zeros(1,sum(npd));
    DT_s=diff(t_segment);
    for i=1:length(npd)
        DT_seg(1+sum(npd(1:(i-1))):sum(npd(1:i)))=repmat(DT_s(i),1,npd(i));
    end
    
    DT_seg=DT_seg';
    
    for i=1:nx
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            cl=[cl,AxLGR*x(:,i)-0.5*DT_seg.*data.data.xrl(i)];
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            ce=[ce,AxLGR*x(:,i)-0.5*DT_seg.*data.data.xrl(i)];
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            cu=[cu,AxLGR*x(:,i)-0.5*DT_seg.*data.data.xru(i)];
        end
    end
    for i=1:nu
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            cl=[cl,AxLGR*[u(:,i);u(end,i)]-0.5*DT_seg.*data.data.url(i)];
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            ce=[ce,AxLGR*[u(:,i);u(end,i)]-0.5*DT_seg.*data.data.url(i)];
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            cu=[cu,AxLGR*[u(:,i);u(end,i)]-0.5*DT_seg.*data.data.uru(i)];
        end
    end
    
    cl=cl(:)';
    cu=cu(:)';
    ce=ce(:)';

else
    DT=diff(t);
    N=size(u,1);
    Ax=diag(-ones(M,1))+diag(ones(M-1,1),1);
    Ax(end,:)=[];
    Au=diag(-ones(N,1))+diag(ones(N-1,1),1);
    Au(end,:)=[];
    for i=1:nx
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            cl=[cl,Ax*x(:,i)-DT.*data.data.xrl(i)];
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            ce=[ce,Ax*x(:,i)-DT.*data.data.xrl(i)];
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            cu=[cu,Ax*x(:,i)-DT.*data.data.xru(i)];
        end
    end
    for i=1:nu
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            cl=[cl,Au*u(:,i)-DT.*data.data.url(i)];
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            ce=[ce,Au*u(:,i)-DT.*data.data.url(i)];
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            cu=[cu,Au*u(:,i)-DT.*data.data.uru(i)];
        end
    end
    cl=reshape(cl',size(cl,1)*size(cl,2),1)';
    cu=reshape(cu',size(cu,1)*size(cu,2),1)';
    ce=reshape(ce',size(ce,1)*size(ce,2),1)';

end

cr=[cl, cu, ce];



end

