function [ solution ] = getVariableRate( solution,data )
%getVariableRate - Obtain the variable rate of change from the discretization mesh.
%
% Syntax:  [ solution ] = getVariableRate( solution,data )
%
% Inputs:
%    solutions - Structure containing the solution
%    data - structure containing matrices to format data 
% 
% Output:
%    solutions - Structure containing the solution with rate data included
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



if (strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))

    t=solution.T;
    x=solution.X;
    u=solution.U;
    
else
    [nt,~,n,m,~,~,M,N,ns]=deal(data.sizes{1:9});
    

    z=solution.z_orgscale;
    Vx=data.map.Vx;
    Vu=data.map.Vu;

    % Generate time vector
    if nt==1
        t0=data.t0;
        tf=solution.tf;
    elseif nt==2
        t0=solution.t0;
        tf=solution.tf;
    else
        t0=data.t0;
        tf=data.tf(1);
    end

    t=(tf-t0)*[0;cumsum(data.tau)]*data.Nm/ns+t0;


    % Extract states and controls

    x=reshape(Vx*z,n,M)';
    u=reshape(Vu*z,m,N)';
end


nx=size(x,2);
nu=size(u,2);
M=size(x,1);


if strcmp(data.options.discretization,'hermite')

    DT=diff(t(1:2:end));

    AxHS1=diag(-3*ones(M,1))+diag(4*ones(M-1,1),1)+diag(-ones(M-2,1),2);
    AxHS1([2:2:end,end],:)=[];
    AuHS1=diag(-3*ones(N,1))+diag(4*ones(N-1,1),1)+diag(-ones(N-2,1),2);
    AuHS1([2:2:end,end],:)=[];
    AxHS3=diag(-ones(M,1))+diag(ones(M-2,1),2);
    AxHS3([2:2:end,end],:)=[];
    AuHS3=diag(-ones(N,1))+diag(ones(N-2,1),2);
    AuHS3([2:2:end,end],:)=[];
    
    dX=zeros(length(DT),nx);dU=zeros(length(DT),nu);
    dX_mid=zeros(length(DT),nx);dU_mid=zeros(length(DT),nu);
    for i=1:nx
        dX(:,i)=AxHS1*x(:,i)./DT;
        dX_mid(:,i)=AxHS3*x(:,i)./DT;
    end
    for i=1:nu
        dU(:,i)=AuHS1*u(:,i)./DT;
        dU_mid(:,i)=AuHS3*u(:,i)./DT;
    end
    dX_org=reshape([dX dX_mid]',size(dX,2),size(dX,1)*2)';
    dU_org=reshape([dU dU_mid]',size(dU,2),size(dU,1)*2)';
    solution.org.dX=dX_org;
    solution.org.dU=dU_org;

elseif (strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))
    
    t0=t(1);
    tf=t(end);
    
%     if (strcmp(data.options.resultRep,'res_min'))
%         AxLGR=data.minres_auxdata.D_mat_quad;
%         npd=data.minres_auxdata.npd_quad+1;
%     else
        AxLGR=data.map.D_structure;
        npd=data.sizes{10};
%     end
    t_segment=(tf-t0)/2*data.tau_segment+(tf+t0)/2;
    DT_seg=zeros(1,sum(npd));
    DT_s=diff(t_segment);
    for i=1:length(npd)
        DT_seg(1+sum(npd(1:(i-1))):sum(npd(1:i)))=repmat(DT_s(i),1,npd(i));
    end
    
    DT_seg=DT_seg';
    
    dX=zeros(length(DT_seg),nx);dU=zeros(length(DT_seg),nu);
    for i=1:nx
        dX(:,i)=AxLGR*x(:,i)./DT_seg./0.5;
    end
    for i=1:nu
        dU(:,i)=AxLGR*[u(:,i);u(end,i)]./DT_seg./0.5;
    end
    

    
else
    DT=diff(t);
    N=size(u,1);
    Ax=diag(-ones(M,1))+diag(ones(M-1,1),1);
    Ax(end,:)=[];
    Au=diag(-ones(N,1))+diag(ones(N-1,1),1);
    Au(end,:)=[];
    
    dX=zeros(length(DT),nx);dU=zeros(length(DT),nu);
    for i=1:nx
            dX(:,i)=Ax*x(:,i)./DT;
    end
    for i=1:nu
            dU(:,i)=Au*u(:,i)./DT;
    end
end


solution.dX=dX;
solution.dU=dU;


end

