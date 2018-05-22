function data = transcriptionMatrixRC( data )
%TRANSCRIPTIONMATRIX -  Precalculate the jacobian calculation of rate
%constraints
%
% Syntax:  data=transcriptionMatrix(data)
%
% Inputs: described in main file
%
% Outputs:
%    data - Structure constaining matrices and vectors
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

if strcmp(data.options.transcription,'globalLGR') || strcmp(data.options.transcription,'hpLGR')
    [nt,np,n,m,~,~,M,N,~,~,~,~,~,~,~,~]=deal(data.sizes{:});
    nz=nt+np+(M+1)*n+M*m;  
else
    [nt,np,n,m,~,~,M,N,~,~,~,~]=deal(data.sizes{:});
    nz=nt+np+M*n+N*m;
end


if strcmp(data.options.transcription,'hermite')
    
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
    
    AxHS1=repelem(AxHS1,1,n);
    AxHS2=repelem(AxHS2,1,n);
    AxHS3=repelem(AxHS3,1,n);
    AxHS4=repelem(AxHS4,1,n);
    AuHS1=repelem(AuHS1,1,m);
    AuHS2=repelem(AuHS2,1,m);
    AuHS3=repelem(AuHS3,1,m);
    AuHS4=repelem(AuHS4,1,m);
   
    
    Acl=[];Ae=[];Acu=[];

    for i=1:n
        idx=linspace(i,i+(M-1)*n,M);
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            AxHScl1=zeros(size(AxHS1));AxHScl1(:,idx)=AxHS1(:,idx);
            AxHScl2=zeros(size(AxHS2));AxHScl2(:,idx)=AxHS2(:,idx);
            AxHScl3=zeros(size(AxHS3));AxHScl3(:,idx)=AxHS3(:,idx);
            Acl=[Acl AxHScl1*data.map.Vx AxHScl2*data.map.Vx AxHScl3*data.map.Vx];
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            AxHSe1=zeros(size(AxHS3));AxHSe1(:,idx)=AxHS3(:,idx);
            AxHSe2=zeros(size(AxHS4));AxHSe2(:,idx)=AxHS4(:,idx);
            Ae=[Ae AxHSe1*data.map.Vx AxHSe2*data.map.Vx];
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            AxHScu1=zeros(size(AxHS1));AxHScu1(:,idx)=AxHS1(:,idx);
            AxHScu2=zeros(size(AxHS2));AxHScu2(:,idx)=AxHS2(:,idx);
            AxHScu3=zeros(size(AxHS3));AxHScu3(:,idx)=AxHS3(:,idx);
            Acu=[Acu AxHScu1*data.map.Vx AxHScu2*data.map.Vx AxHScu3*data.map.Vx];
        end
    end
    for i=1:m
        idx=linspace(i,i+(M-1)*m,M);
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            AuHScl1=zeros(size(AuHS1));AuHScl1(:,idx)=AuHS1(:,idx);
            AuHScl2=zeros(size(AuHS2));AuHScl2(:,idx)=AuHS2(:,idx);
            Acl=[Acl AuHScl1*data.map.Vu AuHScl2*data.map.Vu];
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            AuHSe1=zeros(size(AuHS3));AuHSe1(:,idx)=AuHS3(:,idx);
            AuHSe2=zeros(size(AuHS4));AuHSe2(:,idx)=AuHS4(:,idx);
            Ae=[Aue AuHSe1*data.map.Vu AuHSe2*data.map.Vu];
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            AuHScu1=zeros(size(AuHS1));AuHScu1(:,idx)=AuHS1(:,idx);
            AuHScu2=zeros(size(AuHS2));AuHScu2(:,idx)=AuHS2(:,idx);
            Acu=[Acu AuHScu1*data.map.Vu AuHScu2*data.map.Vu];
        end
    end

    Acl=Acl';
    Acl=reshape(Acl(:),nz,length(Acl(:))/nz)';
    Ae=Ae';
    Ae=reshape(Ae(:),nz,length(Ae(:))/nz)';
    Acu=Acu';
    Acu=reshape(Acu(:),nz,length(Acu(:))/nz)';
    
    data.map.Acl=sparse(Acl);
    data.map.Ae=sparse(Ae);
    data.map.Acu=sparse(Acu);
elseif strcmp(data.options.transcription,'globalLGR') || strcmp(data.options.transcription,'hpLGR')
    AxLGR=data.map.D_structure;
    AuLGR=AxLGR(:,1:end-1);
    AuLGR(:,end)=AuLGR(:,end)+AxLGR(:,end);
   
    Acl=[];Ae=[];Acu=[];

    for i=1:n
        idx=linspace(i,i+(M-1)*n,M);
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            Acl=[Acl;[zeros(size(AxLGR,1),size(AxLGR,2)*(i-1)) AxLGR zeros(size(AxLGR,1),size(AxLGR,2)*(n-i))]*data.map.Vx];
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            Ae=[Ae;[zeros(size(AxLGR,1),size(AxLGR,2)*(i-1)) AxLGR zeros(size(AxLGR,1),size(AxLGR,2)*(n-i))]*data.map.Vx];
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            Acu=[Acu;[zeros(size(AxLGR,1),size(AxLGR,2)*(i-1)) AxLGR zeros(size(AxLGR,1),size(AxLGR,2)*(n-i))]*data.map.Vx];
        end
    end
    for i=1:m
        idx=linspace(i,i+(M-1)*m,M);
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            Acl=[Acl;[zeros(size(AuLGR,1),size(AuLGR,2)*(i-1)) AuLGR zeros(size(AuLGR,1),size(AuLGR,2)*(m-i))]*data.map.Vu];
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            Ae=[Ae;[zeros(size(AuLGR,1),size(AuLGR,2)*(i-1)) AuLGR zeros(size(AuLGR,1),size(AuLGR,2)*(m-i))]*data.map.Vu];
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            Acu=[Acu;[zeros(size(AuLGR,1),size(AuLGR,2)*(i-1)) AuLGR zeros(size(AuLGR,1),size(AuLGR,2)*(m-i))]*data.map.Vu];
        end
    end

    data.map.Acl=sparse(Acl);
    data.map.Ae=sparse(Ae);
    data.map.Acu=sparse(Acu);
else 
    Ax=diag(-ones(M,1))+diag(ones(M-1,1),1);
    Ax(end,:)=[];
    Au=diag(-ones(N,1))+diag(ones(N-1,1),1);
    Au(end,:)=[];
   
    Ax=repelem(Ax,1,n);
    Au=repelem(Au,1,m);
    
    Acl=[];Ae=[];Acu=[];

    for i=1:n
        idx=linspace(i,i+(M-1)*n,M);
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            Axcl1=zeros(size(Ax));Axcl1(:,idx)=Ax(:,idx);
            Acl=[Acl Axcl1*data.map.Vx];
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            Axe1=zeros(size(Ax));Axe1(:,idx)=Ax(:,idx);
            Ae=[Ae Axe1*data.map.Vx];
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            Axcu1=zeros(size(Ax));Axcu1(:,idx)=Ax(:,idx);
            Acu=[Acu Axcu1*data.map.Vx];
        end
    end
    for i=1:m
        idx=linspace(i,i+(M-1)*m,M);
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            Aucl1=zeros(size(Au));Aucl1(:,idx)=Au(:,idx);
            Acl=[Acl Aucl1*data.map.Vu];
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            Aue1=zeros(size(Au));Aue1(:,idx)=Au(:,idx);
            Ae=[Aue Aue1*data.map.Vu];
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            Aucu1=zeros(size(Au));Aucu1(:,idx)=Au(:,idx);
            Acu=[Acu Aucu1*data.map.Vu];
        end
    end

    Acl=Acl';
    Acl=reshape(Acl(:),nz,length(Acl(:))/nz)';
    Ae=Ae';
    Ae=reshape(Ae(:),nz,length(Ae(:))/nz)';
    Acu=Acu';
    Acu=reshape(Acu(:),nz,length(Acu(:))/nz)';
    
    data.map.Acl=sparse(Acl);
    data.map.Ae=sparse(Ae);
    data.map.Acu=sparse(Acu);
end

end

