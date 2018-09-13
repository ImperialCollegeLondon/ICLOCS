function solution=directCollocation(required,z,data)
%DIRECTCOLLOCATION - Generate the cost, constraint and gradient
%information for the direct transcription formulation
%
% Syntax:  solution=directCollocation(required,z,data)
%
% Inputs:
%    required - Flag that determines what to compute for current z
%    z - Current NLP variable
%    data - Structure of data required to compute
%
% Outputs:
%    solution - Data structure containing the solution
%
% Other m-files required: jacobianFD, hessian FD (optional)
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

%------------- BEGIN CODE --------------

%%
global sol;



% Define some useful variables
[nt,np,n,m,ng,~,M,N,ns,~,~,~,~]=deal(data.sizes{:});


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
%     if isfield(vdat.mode, 'np') && vdat.mode.np~=0
%         npactual=np-vdat.mode.np*ng;
%         P=[reshape(repmat(z(nt+1:nt+npactual),M,1),npactual,M)',reshape(z(nt+npactual+1:nt+np),M,ng)];
%     else
        P=reshape(repmat(z(nt+1:nt+np),M,1),np,M)';
%     end
else
    P=spalloc(M,0,1);
end

% Construct time vector
if nt
    tf=z(1);
else
    tf=data.tf; 
end

t0=data.t0;
k0=data.k0;

if data.options.scaling
    t0=scale_variables_back( t0, data.data.Tscale, data.data.Tshift );
    k0=scale_variables_back( k0, data.data.Tscale, data.data.Tshift );
    tf=scale_variables_back( tf, data.data.Tscale, data.data.Tshift );
    data.t0=t0;
    data.k0=k0;
end


% if strcmp(data.options.transcription,'discrete'); tf=1; t0=0; end
T=[0;cumsum(data.tau)]*data.Nm/ns;
t=(tf-t0)*T+k0;

% Extract x0,u0,xf,uf,p
p=z(nt+1:nt+np);
x0=z(nt+np+1:nt+np+n);
u0=z(nt+np+(M)/N*n+1:nt+np+(M)/N*n+m);
xf=z(end-m-n+1:end-m);
uf=z(end-m+1:end);

if strcmp(data.options.derivatives,'adigator')
    [ fx, fxx ] = getAdigator4ICLOCS( X, U, t, p, data );
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



%--------------------------------------------------------------------------
% Return the relevant data
%--------------------------------------------------------------------------


% required
switch required
    
    case{'cost'}
        snm=ones(M,1); 
        sol.cost= mp.w'*((tf-t0)*L(X,Xr,U,Ur,P,t,vdat).*snm)+E(x0,xf,u0,uf,p,t0,tf,vdat);
        solution=sol.cost;
   
    
    case{'gradCost'}
        [sol.gradCost,sol.JL]=gradientCost(L,X,Xr,U,Ur,P,T,E,x0,xf,u0,uf,p,tf,data);
        solution=sol.gradCost;
      
    case{'const'}
        X_vect=reshape(X',n*M,1);
        g_vect=reshape(g(X,U,P,t,vdat)',M*ng,1);
        sol.const=[(x0-data.x0t)*data.cx0;mp.A*X_vect+mp.B*reshape((tf-t0)*f(X,U,P,t,vdat)',M*n,1);
                  g_vect(data.gAllidx);
                  avrc(X,U,P,t,data)';
                  b(x0,xf,u0,uf,p,t0,tf,vdat)];
        solution=sol.const;
        
    case{'jacConst'}
        if strcmp(data.options.derivatives,'numeric')
            sol.jacConst=jacobianFD(f,g,avrc,X,U,P,T,b,x0,xf,u0,uf,p,t0,tf,data);
        end
        if strcmp(data.options.derivatives,'adigator')
            sol.jacConst=jacobianAdigator(f,g,avrc,X,U,P,T,b,x0,xf,u0,uf,p,t0,tf,fx,data);
        end
        if strcmp(data.options.derivatives,'analytic')
            [df,dg,db]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
            [sol.jacConst,sol.Jf]=jacConstz(df,dg,g,f,avrc,X,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data);
        end
          
        if strcmp(data.options.NLPsolver,'worhp')
            solution=sol.jacConst(sub2ind(size(sol.jacConst),data.jacStruct_GRow,data.jacStruct_GCol));
        else
            solution=sol.jacConst;
        end
        
    case{'hessian'}
      if strcmp(data.options.derivatives,'analytic')
            [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,sol.Jf,sol.JL,X,U,P,T,E,b,x0,...
                                                    xf,u0,uf,p,tf,data);                      
            hessc=data.sigma*(Lzz+Ezz)+fzz+gzz+bzz;
            hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
            sol.hessian=tril(hessc);        
      elseif strcmp(data.options.derivatives,'adigator')
            sol.hessian=hessianCDAdigator(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,fxx,data); 
      else
          if strcmp(data.options.hessianFD,'central')
            sol.hessian=hessianCD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data); 
          else
            sol.hessian=hessianFD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data);  
          end
      end

      if strcmp(data.options.NLPsolver,'builtinSQP')
        sol.hessian=sol.hessian+sol.hessian'-diag([diag(sol.hessian)]);
      end
      if strcmp(data.options.NLPsolver,'worhp')
          solution=sol.hessian(sub2ind(size(sol.hessian),data.jacStruct_HMRow,data.jacStruct_HMCol));
      else
          solution=sol.hessian;
      end
end


