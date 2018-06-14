function data=transcriptionMatrix(method,nx,nu,tau,data)
%TRANSCRIPTIONMATRIX -  Format matrices for transcription method
%c(z)=A.Vx.z+B.F(z) and generate the quadrature vector w(tau)
%
% Syntax:  data=transcriptionMatrix(method,nx,nu,tau,data)
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

%------------- BEGIN CODE --------------

% Gets sizes
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{1:9});
nz=nt+np+nx+nu;

if strcmp(method,'globalLGR') || strcmp(method,'hpLGR')

    % Create matrix to extract x vector from z 
    Vx = sparse(nx,nz);
    Vx(1:nx,1:nx)=speye(nx);
    xV=Vx';

    Vx_LGRpts = zeros(nx-n,nz);
    for i=1:n
        Vx_LGRpts(((i-1)*M+1):((i-1)*M+M),((i-1)*M+i):((M+1)*i-1))=speye(M);
    end
    Vx_LGRpts=sparse(Vx_LGRpts);
    xV_LGRpts=Vx_LGRpts';

    % Create matrix to extract u vector from z
    Vu = sparse(nu,nz);
    Vu(1:nu,nx+1:nx+nu)=speye(nu);   
    uV=Vu';

    % Create matrix to extract t vector from z
    Vt = sparse(nt,nz);
    Vt(1:nt,nz-nt+1:nz)=speye(nt);   
    tV=Vt';


    % Create matrices to extract x and u from (x,u) 
    vx=Vx;vu=Vu;
    if nt; vx(:,1)=[];vu(:,1)=[];           end
    if np; vx(:,1:np)=[];vu(:,1:np)=[];     end

else % h methods

    % Preallocate some memory
    A=spalloc(nx-n,nx,3*M*n);B=spalloc(nx-n,nx,3*M*n);w=zeros(M,1);
    Mi=M/N;

    % Create matrix to extract x vector from z 
    Vx=[zeros(nx,nt+np) kron(speye(N),[speye(n*Mi) zeros(n*Mi,m)])];
    xV=Vx';

    % Create matrix to extract u vector from z
    Vu=[zeros(nu,nt+np) kron(speye(N),[zeros(m,n*Mi) speye(m)])];   
    uV=Vu';

    % Create matrices to extract x and u from (x,u) 
    vx=Vx;vu=Vu;
    if nt; vx(:,1)=[];vu(:,1)=[];           end
    if np; vx(:,1:np)=[];vu(:,1:np)=[];     end

end


% Generate A,B matrices and quadrature vector w
switch method

    case{'discrete'}

        disp('Formatting matrices for the discrete-time system')

        A=spdiags(ones(M*n,1),n,(M-1)*n,M*n);

        B=-spdiags(ones(n*M,1),0,n*(M-1),n*M);

        w=[ones(M-1,1);0];
        W=sparse(1:M,1:M,w);

    case{'euler'}

        disp('Formatting matrices for the euler approximation')
        d1=spdiags(-ones(M*n,1),0,(M-1)*n,M*n);
        d2=spdiags(ones(M*n,1),n,(M-1)*n,M*n);
        A=d1+d2;

       % h=repmat(tau,n,1);
        
        h=kron(tau,ones(n,1));   %PF
        B=-spdiags(h,0,(M-1)*n,n*M);
        
        w=[tau;0];
        W=sparse(1:M,1:M,w);

    case{'trapezoidal'}

        disp('Formatting matrices for the trapezoidal approximation')
        d1=spdiags(-ones(M*n,1),0,(M-1)*n,M*n);
        d2=spdiags(ones(M*n,1),n,(M-1)*n,M*n);
        
        A=d1+d2;

     

       h1=-0.5*kron(tau,ones(n,1));
       B=spdiags(h1,0,(M-1)*n,n*M)+spdiags(h1,n,(M-1)*n,n*M);


        w=0.5*[tau(1);tau(1:M-2)+tau(2:M-1);tau(M-1)];   % PF
        W=sparse(1:M,1:M,w);    

     
        
        
      case{'hermite'}

        disp('Formatting matrices for the hermite-simpson approximation')
        A=spdiags(-0.5*kron(ones(M,1),[ones(n,1);zeros(n,1)]),2*n,A);
        A=spdiags(ones((2*M)*n,1),n,A);
        A=spdiags(-0.5*kron(ones(M,1),[ones(n,1);zeros(n,1)]),0,A);
        A=spdiags(-kron(ones(M-1,1),[zeros(n,1);ones(n,1)]),-n,A);
        
      
        t0 =-tau(:,ones(1,n)).';ta = t0(:)/8;tb=t0(:)/6;tc=2*t0(:)/3;
        B=spdiags(ta.*repmat([ones(n,1);zeros(n,1)],(M-1)/2,1),0,B);
        B=spdiags(-ta.*repmat([ones(n,1);zeros(n,1)],(M-1)/2,1),2*n,B);
        B=spdiags(tb.*repmat([zeros(n,1);ones(n,1)],(M-1)/2,1),-n,B);
        B=spdiags(tb.*repmat([zeros(n,1);ones(n,1)],(M-1)/2,1),n,B);
        B=B+spdiags(tc.*repmat([zeros(n,1);ones(n,1)],(M-1)/2,1),0,nx-n,nx);

        vtau=tau(1:2:M-1);
        w0(1:2:M-2,1)=4*vtau;
        w0(2:2:M-2,1)=vtau(1:end-1)+vtau(2:end);
        
        w=[tau(1);w0;tau(end)]/6;
        W=sparse(1:length(w),1:length(w),w);

      case{'globalLGR','hpLGR'}

        disp('Formatting matrices for the LGR quadrature')
        
    otherwise
        disp('Method unknown. Please check spelling')
end



% Store the matrices and return to formatNLP.m

map.Vx=Vx;map.xV=xV;map.Vu=Vu;map.uV=uV;


if strcmp(method,'globalLGR') || strcmp(method,'hpLGR')
    map.Vx_LGRpts=Vx_LGRpts;map.xV_LGRpts=xV_LGRpts;map.Vt=Vt;map.tV=tV;
else
    map.A=A;map.B=B;map.w=w;map.W=W;
end

data.map=map;