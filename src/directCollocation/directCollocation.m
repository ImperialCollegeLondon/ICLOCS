function [solution,XU0f]=directCollocation(required,z,data,phaseNo)
%DIRECTCOLLOCATION - Generate the cost, constraint and gradient
%information for the direct transcription formulation
%
% Syntax:  solution=directCollocation(required,z,data)
%
% Inputs:
%    required - Flag that determines what to compute for current z
%    z - Current NLP variable
%    data - Structure of data required to compute
%    phaseNo - phase number
%
% Outputs:
%    solution - Data structure containing the solution
%    XU0f - data for linkage constraints
%
% Other m-files required: jacobianFD, hessian FD (optional)
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

%------------- BEGIN CODE --------------

%%
global sol;



% Define some useful variables
[nt,np,n,m,ng,~,M,N,ns,~,~,~,~]=deal(data.sizes{1:13});


% Get function definitions
[L,E,f,g,avrc,b]=deal(data.functions{:});
vdat=data.data;

% Get matrices
mp=data.map;

%--------------------------------------------------------------------------
% Extract and format vectors from NLP variable
%--------------------------------------------------------------------------

% Extract states and inputs from z and reshape for function evaluations
X=reshape(mp.Vx*z,n,M)';
U=reshape(mp.Vu*z,m,N)';


%%
% Extract design parameters if specified and convert to cells
if np
    P=reshape(repmat(z(nt+1:nt+np),M,1),np,M)';
else
    P=spalloc(M,0,1);
end

% Construct time vector
if nt==1
    tf=z(1);
    t0=data.t0;
elseif nt==2
    t0=z(1);
    tf=z(2);
else
    tf=data.tf; 
    t0=data.t0;
end



% if strcmp(data.options.transcription,'discrete'); tf=1; t0=0; end
tau=[0;cumsum(data.tau)]*data.Nm/ns;
t=(tf-t0)*tau+t0;
% t=(tf-t0)*tau+k0;

% Extract x0,u0,xf,uf,p
p=z(nt+1:nt+np);
x0=z(nt+np+1:nt+np+n);
u0=z(nt+np+(M)/N*n+1:nt+np+(M)/N*n+m);
xf=z(end-m-n+1:end-m);
uf=z(end-m+1:end);

if strcmp(data.options.derivatives,'adigator') 
    if strcmp(data.options.discretization,'discrete') 
        [const_vec_Adigator] = getAdigator4ICLOCS( X, U, linspace(0,1,length(t)), p, data );
    else
        [const_vec_Adigator] = getAdigator4ICLOCS( X, U, t, p, data );
    end
end

% Format reference inputs and states if applicable
if ~isempty(data.references.xr)
    Xr=data.references.xr;
else
    Xr=[];
end

if ~isempty(data.references.ur)
    Ur=data.references.ur;
else
    Ur=[];
end

% Prepare the data for linkage constraints
[XU0f] = prepareForLinkageConstraintInfo(x0,xf,u0,uf,vdat);

%--------------------------------------------------------------------------
% Return the relevant data
%--------------------------------------------------------------------------


% required
switch required
    
    case{'cost'}
        snm=ones(M,1); 
        sol{phaseNo}.cost= mp.w'*((tf-t0)*L(X,Xr,U,Ur,P,t,vdat).*snm)+E(x0,xf,u0,uf,p,t0,tf,vdat);
        solution=sol{phaseNo}.cost;
   
    
    case{'gradCost'}
        [sol{phaseNo}.gradCost,sol{phaseNo}.JL]=gradientCost(L,X,Xr,U,Ur,P,tau,E,x0,xf,u0,uf,p,t0,tf,data);
        solution=sol{phaseNo}.gradCost;
      
    case{'const'}
        if strcmp(data.options.derivatives,'adigator')
            sol{phaseNo}.const=[const_vec_Adigator.f;
                      avrc(X,U,P,t,data)';
                      b(x0,xf,u0,uf,p,t0,tf,vdat)];
        else
            X_vect=reshape(X',n*M,1);
            g_vect=reshape(g(X,U,P,t,vdat)',M*ng,1);
            sol{phaseNo}.const=[(x0-data.x0t)*data.cx0;mp.A*X_vect+mp.B*reshape((tf-t0)*f(X,U,P,t,vdat)',M*n,1);
                      g_vect(data.gAllidx);
                      avrc(X,U,P,t,data)';
                      b(x0,xf,u0,uf,p,t0,tf,vdat)];
        end
        solution=sol{phaseNo}.const;
        
    case{'jacConst'}
        if strcmp(data.options.derivatives,'numeric')
            sol{phaseNo}.jacConst=jacobianFD(f,g,avrc,X,U,P,tau,b,x0,xf,u0,uf,p,t0,tf,data);
        end
        if strcmp(data.options.derivatives,'adigator')
            sol{phaseNo}.jacConst=jacobianAdigator(f,g,avrc,X,U,P,tau,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data);
        end
        if strcmp(data.options.derivatives,'analytic')
            jacConst=data.analyticDeriv.jacConst;
            [df,dg,db]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
            [sol{phaseNo}.jacConst,sol{phaseNo}.Jf]=jacConstz(df,dg,g,f,avrc,X,U,P,tau,db,b,x0,xf,u0,uf,p,t0,tf,data);
        end
        if strcmp(data.options.NLPsolver,'worhp')
            solution=sol{phaseNo}.jacConst(sub2ind(size(sol{phaseNo}.jacConst),data.jacStruct_GRow,data.jacStruct_GCol));
        else
            solution=sol{phaseNo}.jacConst;
        end
        
    case{'hessian'}
      if strcmp(data.options.derivatives,'analytic')
            [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,sol{phaseNo}.Jf,sol{phaseNo}.JL,X,U,P,tau,E,b,x0,...
                                                    xf,u0,uf,p,t0,tf,data);                      
            hessc=data.sigma*(Lzz+Ezz)+fzz+gzz+bzz;
            hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
            sol{phaseNo}.hessian=tril(hessc);        
      elseif strcmp(data.options.derivatives,'adigator')
            sol{phaseNo}.hessian=hessianCDAdigator(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data); 
      else
            sol{phaseNo}.hessian=hessianCD(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data); 
      end

      if strcmp(data.options.NLPsolver,'builtinSQP')
        sol{phaseNo}.hessian=sol{phaseNo}.hessian+sol{phaseNo}.hessian'-diag([diag(sol{phaseNo}.hessian)]);
      end
      if strcmp(data.options.NLPsolver,'worhp')
          solution=sol{phaseNo}.hessian(sub2ind(size(sol{phaseNo}.hessian),data.jacStruct_HMRow,data.jacStruct_HMCol));
      else
          solution=sol{phaseNo}.hessian;
      end
end


