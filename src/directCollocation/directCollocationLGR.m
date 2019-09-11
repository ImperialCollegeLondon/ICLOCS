function [solution,XU0f]=directCollocationLGR(required,z,data,phaseNo)
%DIRECTCOLLOCATION - Generate the cost, constraint and gradient
%information for the LQR direct transcription formulation
%
% Syntax:  [solution,XU0f]=directCollocationLGR(required,z,data,phaseNo)
%
% Inputs:
%    required - Flag that determines what to compute for current z
%    z - Current NLP variable
%    data - Structure of data required to compute
%    phaseNo - phase number
% 
% Outputs:
%    solution - Data structure containing the solution
%
% Other m-files required: jacobianFD_LGR, hessianCD_LGR (optional)
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
global sol ro_time;

% Define some useful variables
[nt,np,n,m,ng,~,M,~,~,npd,~,npduidx,nps,~,~,~,~]=deal(data.sizes{1:17});

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
    t0=z(end-nt+1);
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
    data.tau_segment=tau_segment';
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
uf=U(end,:)'; %End point input obtained via extrapolation

if strcmp(data.options.derivatives,'adigator')
    [const_vec_Adigator] = getAdigator4ICLOCS( X_Np1, U, [t;tf], p, data );
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
t_segment_end=data.t_segment_mat_m*t_segment;
k0=data.t_segment_mat_p*t_segment;
data.t_segment_end=t_segment_end;
data.diag_t_segment_end=sparse(diag(t_segment_end));
data.k0=k0;
        
D_structure=data.map.D_structure;
w=data.map.w;

% Prepare the data for linkage constraints
[XU0f] = prepareForLinkageConstraintInfo(x0,xf,u0,uf,vdat);
%%
%--------------------------------------------------------------------------
% Return the relevant data
%--------------------------------------------------------------------------



switch required

    case{'cost'} 
        sol{phaseNo}.cost= w'*(t_segment_end.*L(X,Xr,U,Ur,P,t,vdat))+E(x0,xf,u0,uf,p,t0,tf,vdat);
        solution=sol{phaseNo}.cost;
        
    case{'gradCost'}
        [sol{phaseNo}.gradCost,sol{phaseNo}.JL]=gradientCostLGR(L,X,Xr,U,Ur,P,tau_inc,t,E,x0,xf,u0,uf,p,t0,tf,data);
        if data.options.reorderLGR
            solution=sol{phaseNo}.gradCost(data.reorder.z_idx);
        else
            solution=sol{phaseNo}.gradCost;
        end
        
    case{'const'} 
        if strcmp(data.options.derivatives,'adigator')
            sol{phaseNo}.const=[const_vec_Adigator.f;
                       avrc(X_Np1,U,P,[t;tf],data)';
                       b(x0,xf,u0,uf,p,t0,tf,vdat,data.options,t_segment)];
                  
        else
            g_vect=reshape(g(X,U,P,t,vdat),M*ng,1);
            sol{phaseNo}.const=[reshape(D_structure*X_Np1-diag(t_segment_end)*f(X,U,P,t,vdat),M*n,1);
               g_vect(data.gAllidx);
               avrc(X_Np1,U,P,[t;tf],data)';
               b(x0,xf,u0,uf,p,t0,tf,vdat,data.options,t_segment)];
        end

        if data.options.reorderLGR
            solution=sol{phaseNo}.const(data.reorder.vert_idx);
        else
            solution=sol{phaseNo}.const;
        end
        
    case{'jacConst'}
        if strcmp(data.options.derivatives,'numeric')
            sol{phaseNo}.jacConst=jacobianFD_LGR(f,g,avrc,X_Np1,U,P,tau_inc,b,x0,xf,u0,uf,p,t0,tf,data);
        elseif strcmp(data.options.derivatives,'adigator')
            sol{phaseNo}.jacConst=jacobianFDAdigator_LGR(f,g,avrc,X_Np1,U,P,tau_inc,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data);
        elseif strcmp(data.options.derivatives,'analytic')
            jacConst=data.analyticDeriv.jacConst;
            [df,dg,db]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
            [sol{phaseNo}.jacConst,sol{phaseNo}.Jf]=jacConstzLGR(df,dg,g,f,X,U,P,tau_inc,db,b,x0,xf,u0,uf,p,t0,tf,data);
        end
        if data.options.reorderLGR
            sol{phaseNo}.jacConst=sol{phaseNo}.jacConst(data.reorder.vert_idx,data.reorder.z_idx);
        end
        
        
        if strcmp(data.options.NLPsolver,'worhp')
            solution=sol{phaseNo}.jacConst(sub2ind(size(sol{phaseNo}.jacConst),data.jacStruct_GRow,data.jacStruct_GCol));
        else
            solution=sol{phaseNo}.jacConst;
        end
       
    
    case{'hessian'}
      if strcmp(data.options.derivatives,'analytic')
            [Lzz,Ezz,fzz,gzz,bzz]=hessianAN_LGR(L,f,g,sol{phaseNo}.Jf,sol{phaseNo}.JL,X,U,P,tau_inc,E,b,x0,...
                                                    xf,u0,uf,p,t0,tf,data);
               hessc=data.sigma*(Lzz+Ezz)-fzz+gzz+bzz;
               hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
               sol{phaseNo}.hessian=hessc;     
      elseif strcmp(data.options.derivatives,'adigator')
            sol{phaseNo}.hessian=hessianCDAdigator_LGR(L,f,g,X,U,P,tau_inc,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data); 
      else
            sol{phaseNo}.hessian=hessianCD_LGR(L,f,g,X,U,P,tau_inc,E,b,x0,xf,u0,uf,p,t0,tf,data); 
      end

      if data.options.reorderLGR
          sol{phaseNo}.hessian=sol{phaseNo}.hessian(data.reorder.z_idx,data.reorder.z_idx);
      end    
      
      sol{phaseNo}.hessian=tril(sol{phaseNo}.hessian);
      if strcmp(data.options.NLPsolver,'worhp')
          solution=sol{phaseNo}.hessian(sub2ind(size(sol{phaseNo}.hessian),data.jacStruct_HMRow,data.jacStruct_HMCol));
      else
          solution=sol{phaseNo}.hessian;
      end

end



