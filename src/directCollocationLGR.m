function [solution,data]=directCollocationLGR(required,z,data)
%DIRECTCOLLOCATION - Generate the cost, constraint and gradient
%information for the LQR direct transcription formulation
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
% Other m-files required: jacobianFD_LGR, hessianCD_LGR (optional)
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
global sol ro_time;

% Define some useful variables
[nt,np,n,m,ng,~,M,~,~,npd,~,npduidx,nps,~,~,~]=deal(data.sizes{:});

% Reorder for alternative LGR formulation with timing
tic;
if data.options.reorderLGR
    z=z(data.reorder.z_idx_back);
    if (strcmp(required,'hessian'))
        data.lambda=data.lambda(data.reorder.vert_idx_back);
    end
end
elapsedTime = toc;
ro_time=ro_time+elapsedTime;

% Get function definitions
[L,E,f,g,avrc,b]=deal(data.functions{:});
vdat=data.data;

% Get matrices
mp=data.map;

%--------------------------------------------------------------------------
% Extract and format vectors from NLP variable
%--------------------------------------------------------------------------

% Extract states and inputs from z and reshape for function evaluations
X_Np1=reshape(mp.Vx*z,M+1,n); %States including the end point
X=X_Np1(1:M,:); %States exlduing the end point
U=reshape(mp.Vu*z,M,m); %input exluding the end point

%%
% Extract design parameters if specified and convert to cells
if np
    P=reshape(repmat(z(n*(M+1)+M*m+1:n*(M+1)+M*m+np),M,1),M,np);
else
    P=spalloc(M,0,1);
end

% Construct time vector
if (nt==1) 
    tf=z(1);
    t0=data.t0;
elseif (nt>=2) %LGR
    t0=data.t0;
    tf=z(end);
else
    tf=data.tf; 
    t0=data.t0;
end


if data.options.adaptseg==1 %Adaptive LGR method
    %Update the mesh with information from optimization time variables
    t_segment=z(end-nt+1:end);
    tau_segment = ((1-(-1))*t_segment-(t_segment(1)+t_segment(end))) /(t_segment(end)-t_segment(1));
    tau_inc=[];
    for i=1:length(tau_segment)-1
        tau_seg{i} = ((tau_segment(i+1)+tau_segment(i)) + (tau_segment(i+1)-tau_segment(i)) * [data.map.LGR.points{npduidx(i)}; 1]) / (1-(-1));
        tau_inc=[tau_inc; tau_seg{i}(1:end-1)];
    end
    tau=diff(tau_inc);
    data.tau_inc=tau_inc;
    data.t_segment=t_segment;
    data.tau_segment=tau_segment;
    data.tau=tau;
    t=(tf-t0)/2*tau_inc+(tf+t0)/2;
else
    tau_inc=data.tau_inc;
    t=(tf-t0)/2*tau_inc+(tf+t0)/2;

    %for multi-interval mesh
    if nt>=2
        t_segment=(tf-t0)/2*data.tau_segment'+(tf+t0)/2; %Time at start/end of each segment
        data.t_segment=[t_segment(1); t_segment(end)];
    end
end



% Extract x0,u0,xf,uf,p
p=z(n*(M+1)+M*m+1:n*(M+1)+M*m+np);
x0=X_Np1(1,:)';
u0=U(1,:)';
xf=X_Np1(end,:)';
uf=U(end,:)';
% uf=interp1(tau_inc,U,1,'linear','extrap')'; %End point input obtained via extrapolation

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

% Generate transformation matrices
t_segment_mat_m=zeros(M,nt);
prow=cumsum([1,npd]);
for i=1:nps
    t_segment_mat_m(prow(i):prow(i)+npd(i)-1,i)=-0.5;
    t_segment_mat_m(prow(i):prow(i)+npd(i)-1,i+1)=0.5;
end
t_segment_end=t_segment_mat_m*t_segment;
t_segment_mat_p=zeros(M,nt);
for i=1:nps
    t_segment_mat_p(prow(i):prow(i)+npd(i)-1,i)=0.5;
    t_segment_mat_p(prow(i):prow(i)+npd(i)-1,i+1)=0.5;
end
k0=t_segment_mat_p*t_segment;

data.t_segment_mat_m=t_segment_mat_m;
data.t_segment_mat_p=t_segment_mat_p;
data.t_segment_end=t_segment_end;
data.k0=k0;
        
D_structure=data.map.D_structure;
w=data.map.w;
%%
%--------------------------------------------------------------------------
% Return the relevant data
%--------------------------------------------------------------------------



switch required

    case{'cost'} 
        sol.cost= w'*(t_segment_end.*L(X,Xr,U,Ur,P,t,vdat))+E(x0,xf,u0,uf,p,t0,tf,vdat);
        solution=sol.cost;
        
    case{'gradCost'}
        [sol.gradCost,sol.JL]=gradientCostLGR(L,X,Xr,U,Ur,P,tau_inc,t,E,x0,xf,u0,uf,p,t0,tf,data);
        if data.options.reorderLGR
            tic;
            solution=sol.gradCost(data.reorder.z_idx);
            elapsedTime = toc;
            ro_time=ro_time+elapsedTime;
        else
            solution=sol.gradCost;
        end
        
    case{'const'} 
        sol.const=[reshape(D_structure*X_Np1-diag(t_segment_end)*f(X,U,P,t,vdat),M*n,1);
           reshape(g(X,U,P,t,vdat),M*ng,1);
           avrc(X_Np1,U,P,[t;tf],data)';
           b(x0,xf,u0,uf,p,t0,tf,vdat,data.options,t_segment)];
        if data.options.reorderLGR
            solution=sol.const(data.reorder.vert_idx);
        else
            solution=sol.const;
        end
        
    case{'jacConst'}
        if strcmp(data.options.derivatives,'numeric')
            sol.jacConst=jacobianFD_LGR(f,g,avrc,X_Np1,U,P,tau_inc,b,x0,xf,u0,uf,p,t0,tf,data);
        elseif strcmp(data.options.derivatives,'adigator')
            sol.jacConst=jacobianFDAdigator_LGR(f,g,avrc,X_Np1,U,P,tau_inc,b,x0,xf,u0,uf,p,t0,tf,fx,data);
        elseif strcmp(data.options.derivatives,'analytic')
            [df,dg,db]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
            [sol.jacConst,sol.Jf]=jacConstzLGR(df,dg,g,f,X,U,P,tau_inc,db,b,x0,xf,u0,uf,p,t0,tf,data);
        end
        if data.options.reorderLGR
            sol.jacConst=sol.jacConst(data.reorder.vert_idx,data.reorder.z_idx);
        end
        if strcmp(data.options.NLPsolver,'worhp')
            solution=sol.jacConst(sub2ind(size(sol.jacConst),data.jacStruct_GRow,data.jacStruct_GCol));
        else
            solution=sol.jacConst;
        end
       
    
    case{'hessian'}
      if strcmp(data.options.derivatives,'analytic')
            [Lzz,Ezz,fzz,gzz,bzz]=hessianAN_LGR(L,f,g,sol.Jf,sol.JL,X,U,P,tau_inc,E,b,x0,...
                                                    xf,u0,uf,p,t0,tf,data);
               hessc=data.sigma*(Lzz+Ezz)-fzz+gzz+bzz;
               hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
               sol.hessian=hessc;     
      elseif strcmp(data.options.derivatives,'adigator')
          if strcmp(data.options.hessianFD,'central')
            sol.hessian=hessianCDAdigator_LGR(L,f,g,X,U,P,tau_inc,E,b,x0,xf,u0,uf,p,t0,tf,fxx,data); 
          else
%             sol.hessian=hessianFDAdigator_LGR(L,f,g,X,U,P,tau_inc,E,b,x0,xf,u0,uf,p,t0,tf,fxx,data);  
          end
      else
          if strcmp(data.options.hessianFD,'central')
            sol.hessian=hessianCD_LGR(L,f,g,X,U,P,tau_inc,E,b,x0,xf,u0,uf,p,t0,tf,data); 
          else
%             sol.hessian=hessianFD_LGR(L,f,g,X,U,P,tau_inc,E,b,x0,xf,u0,uf,p,t0,tf,data);  
          end
      end

      if data.options.reorderLGR
          sol.hessian=sol.hessian(data.reorder.z_idx,data.reorder.z_idx);
      end    
      sol.hessian=tril(sol.hessian);
      if strcmp(data.options.NLPsolver,'worhp')
          solution=sol.hessian(sub2ind(size(sol.hessian),data.jacStruct_HMRow,data.jacStruct_HMCol));
      else
          solution=sol.hessian;
      end

end



