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
    [nt,np,nx,nu,~,~,M,N,ns,npd,~,~,~,nrcl,nrcu,nrce,~]=deal(data.sizes{:});
    nz=nt+np+(M+1)*nx+M*nu;  
else
    [nt,np,nx,nu,~,~,M,N,ns,nrcl,nrcu,nrce,~]=deal(data.sizes{:});
    nz=nt+np+M*nx+N*nu;
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
    
    data.data.RC.AxHS1=AxHS1;
    data.data.RC.AxHS2=AxHS2;
    data.data.RC.AxHS3=AxHS3;
    data.data.RC.AxHS4=AxHS4;
    data.data.RC.AuHS1=AuHS1;
    data.data.RC.AuHS2=AuHS2;
    data.data.RC.AuHS3=AuHS3;
    data.data.RC.AuHS4=AuHS4;
    
    AxHS1=repelem(AxHS1,1,nx);
    AxHS2=repelem(AxHS2,1,nx);
    AxHS3=repelem(AxHS3,1,nx);
    AxHS4=repelem(AxHS4,1,nx);
    AuHS1=repelem(AuHS1,1,nu);
    AuHS2=repelem(AuHS2,1,nu);
    AuHS3=repelem(AuHS3,1,nu);
    AuHS4=repelem(AuHS4,1,nu);
   
    Mt=(M-1)/ns;
    Acl=zeros(Mt,nrcl/Mt*nz);
    Ae=zeros(Mt,nrce/Mt*nz);
    Acu=zeros(Mt,nrcu/Mt*nz);
    k=1;l=1;m=1;
    
    for i=1:nx
        idx=linspace(i,i+(M-1)*nx,M);
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            AxHScl1=zeros(size(AxHS1));AxHScl1(:,idx)=AxHS1(:,idx);
            AxHScl2=zeros(size(AxHS2));AxHScl2(:,idx)=AxHS2(:,idx);
            AxHScl3=zeros(size(AxHS3));AxHScl3(:,idx)=AxHS3(:,idx);
            Acl(:,k:k+nz*3-1)=[sparse(AxHScl1)*data.map.Vx sparse(AxHScl2)*data.map.Vx sparse(AxHScl3)*data.map.Vx];
            k=k+nz*3;
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            AxHSe1=zeros(size(AxHS3));AxHSe1(:,idx)=AxHS3(:,idx);
            AxHSe2=zeros(size(AxHS4));AxHSe2(:,idx)=AxHS4(:,idx);
            Ae(:,l:l+nz*2-1)=[sparse(AxHSe1)*data.map.Vx sparse(AxHSe2)*data.map.Vx];
            l=l+nz*2;
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            AxHScu1=zeros(size(AxHS1));AxHScu1(:,idx)=AxHS1(:,idx);
            AxHScu2=zeros(size(AxHS2));AxHScu2(:,idx)=AxHS2(:,idx);
            AxHScu3=zeros(size(AxHS3));AxHScu3(:,idx)=AxHS3(:,idx);
            Acu(:,m:m+nz*3-1)=[sparse(AxHScu1)*data.map.Vx sparse(AxHScu2)*data.map.Vx sparse(AxHScu3)*data.map.Vx];
            m=m+nz*3;
        end
    end
    for i=1:nu
        idx=linspace(i,i+(M-1)*nu,M);
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            AuHScl1=zeros(size(AuHS1));AuHScl1(:,idx)=AuHS1(:,idx);
            AuHScl2=zeros(size(AuHS2));AuHScl2(:,idx)=AuHS2(:,idx);
            Acl(:,k:k+nz*2-1)=[sparse(AuHScl1)*data.map.Vu sparse(AuHScl2)*data.map.Vu];
            k=k+nz*2;
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            AuHSe1=zeros(size(AuHS3));AuHSe1(:,idx)=AuHS3(:,idx);
            AuHSe2=zeros(size(AuHS4));AuHSe2(:,idx)=AuHS4(:,idx);
            Ae(:,l:l+nz*2-1)=[sparse(AuHSe1)*data.map.Vu sparse(AuHSe2)*data.map.Vu];
            l=l+nz*2;
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            AuHScu1=zeros(size(AuHS1));AuHScu1(:,idx)=AuHS1(:,idx);
            AuHScu2=zeros(size(AuHS2));AuHScu2(:,idx)=AuHS2(:,idx);
            Acu(:,m:m+nz*2-1)=[sparse(AuHScu1)*data.map.Vu sparse(AuHScu2)*data.map.Vu];
            m=m+nz*2;
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
   
    Mt=sum(npd);
    Acl=zeros(nrcl,nz);
    Ae=zeros(nrce,nz);
    Acu=zeros(nrcu,nz);
    k=1;l=1;m=1;

    for i=1:nx
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            Acl(k:k+Mt-1,:)=sparse([zeros(size(AxLGR,1),size(AxLGR,2)*(i-1)) AxLGR zeros(size(AxLGR,1),size(AxLGR,2)*(nx-i))])*data.map.Vx;
            k=k+Mt;
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            Ae(l:l+Mt-1,:)=sparse([zeros(size(AxLGR,1),size(AxLGR,2)*(i-1)) AxLGR zeros(size(AxLGR,1),size(AxLGR,2)*(nx-i))])*data.map.Vx;
            l=l+Mt;
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            Acu(m:m+Mt-1,:)=sparse([zeros(size(AxLGR,1),size(AxLGR,2)*(i-1)) AxLGR zeros(size(AxLGR,1),size(AxLGR,2)*(nx-i))])*data.map.Vx;
            m=m+Mt;
        end
    end
    for i=1:nu
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            Acl(k:k+Mt-1,:)=sparse([zeros(size(AuLGR,1),size(AuLGR,2)*(i-1)) AuLGR zeros(size(AuLGR,1),size(AuLGR,2)*(nu-i))])*data.map.Vu;
            k=k+Mt;
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            Ae(l:l+Mt-1,:)=sparse([zeros(size(AuLGR,1),size(AuLGR,2)*(i-1)) AuLGR zeros(size(AuLGR,1),size(AuLGR,2)*(nu-i))])*data.map.Vu;
            l=l+Mt;
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            Acu(m:m+Mt-1,:)=sparse([zeros(size(AuLGR,1),size(AuLGR,2)*(i-1)) AuLGR zeros(size(AuLGR,1),size(AuLGR,2)*(nu-i))])*data.map.Vu;
            m=m+Mt;
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
   
    Ax=repelem(Ax,1,nx);
    Au=repelem(Au,1,nu);
    
    Mt=M-1;
    Acl=zeros(Mt,nrcl/Mt*nz);
    Ae=zeros(Mt,nrce/Mt*nz);
    Acu=zeros(Mt,nrcu/Mt*nz);
    k=1;l=1;m=1;

    for i=1:nx
        idx=linspace(i,i+(M-1)*nx,M);
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)~=data.data.xru(i)
            Axcl1=zeros(size(Ax));Axcl1(:,idx)=Ax(:,idx);
            Acl(:,k:k+nz-1)=sparse(Axcl1)*data.map.Vx;
            k=k+nz;
        end
        if ~isinf(data.data.xrl(i)) && data.data.xrl(i)==data.data.xru(i)
            Axe1=zeros(size(Ax));Axe1(:,idx)=Ax(:,idx);
            Ae(:,l:l+nz-1)=sparse(Axe1)*data.map.Vx;
            l=l+nz;
        end
        if ~isinf(data.data.xru(i)) && data.data.xrl(i)~=data.data.xru(i)
            Axcu1=zeros(size(Ax));Axcu1(:,idx)=Ax(:,idx);
            Acu(:,m:m+nz-1)=sparse(Axcu1)*data.map.Vx;
            m=m+nz;
        end
    end
    for i=1:nu
        idx=linspace(i,i+(M-1)*nu,M);
        if ~isinf(data.data.url(i)) && data.data.url(i)~=data.data.uru(i)
            Aucl1=zeros(size(Au));Aucl1(:,idx)=Au(:,idx);
            Acl(:,k:k+nz-1)=sparse(Aucl1)*data.map.Vu;
            k=k+nz;
        end
        if ~isinf(data.data.url(i)) && data.data.url(i)==data.data.uru(i)
            Aue1=zeros(size(Au));Aue1(:,idx)=Au(:,idx);
            Ae(:,l:l+nz-1)=sparse(Aue1)*data.map.Vu;
            l=l+nz;
        end
        if ~isinf(data.data.uru(i)) && data.data.url(i)~=data.data.uru(i)
            Aucu1=zeros(size(Au));Aucu1(:,idx)=Au(:,idx);
            Acu(:,m:m+nz-1)=sparse(Aucu1)*data.map.Vu;
            m=m+nz;
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

