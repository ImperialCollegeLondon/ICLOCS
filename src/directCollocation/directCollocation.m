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

sol{phaseNo}.z=z;

try

    % Define some useful variables
    [nt,np,n,m,ng,nb,M,N,ns,~,~,~,~]=deal(data.sizes{1:13});
    [ng_eq,ng_neq]=deal(data.sizes{15:16});

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
    if strcmp(data.options.discretization,'discrete') || strcmp(data.options.discretization,'euler')
        U(end,:)=U(end-1,:);
    end
%     if isfield(data,'t_zone_map')
%         U=U(data.t_zone_map,:);
%     end


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
    if strcmp(data.options.discretization,'discrete') 
        tau=[0;data.tau_inc]*(data.Nm-1)/ns;
    else
        tau=[0;data.tau_inc]*data.Nm/ns;
    end
    t=(tf-t0)*tau+t0;
    % t=(tf-t0)*tau+k0;

    % Extract x0,u0,xf,uf,p
    p=z(nt+1:nt+np);
    x0=z(nt+np+1:nt+np+n);
    u0=z(nt+np+(M)/N*n+1:nt+np+(M)/N*n+m);
    xf=z(end-m-n+1:end-m);
    if strcmp(data.options.discretization,'discrete') || strcmp(data.options.discretization,'euler')
        uf=U(end,:)';
    else
        uf=z(end-m+1:end);
    end
    

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

    if  strcmp(data.options.transcription,'direct_collocation_intres_reg')
        intres_reg=1;
        if strcmp(data.options.derivatives,'adigator')
            [ ResNorm_vec, ~] = getAdigator4ICLOCS_minres( X, U, t, P, data.resmin );
        end
    else
        intres_reg=0;
    end

    %--------------------------------------------------------------------------
    % Return the relevant data
    %--------------------------------------------------------------------------


    % required
    switch required

        case{'cost'}
            snm=ones(M,1); 
            if intres_reg
                if strcmp(data.resmin.dataNLP.options.derivatives,'adigator')
                    ResNorm_v=ResNorm_vec.f;
                else
                    [ResNorm] = costResidualMin_ModeMinRes( X,U,P,t,data.resmin);
                    ResNorm_v=sum(ResNorm,2);
                end
                sol{phaseNo}.cost= mp.w'*((tf-t0)*L(X,Xr,U,Ur,P,t,vdat).*snm)+E(x0,xf,u0,uf,p,t0,tf,vdat)+ResNorm_v;
            else
                sol{phaseNo}.cost= mp.w'*((tf-t0)*L(X,Xr,U,Ur,P,t,vdat).*snm)+E(x0,xf,u0,uf,p,t0,tf,vdat);
            end
            solution=sol{phaseNo}.cost;
%             

        case{'gradCost'}
            [sol{phaseNo}.gradCost,sol{phaseNo}.JL]=gradientCost(L,X,Xr,U,Ur,P,tau,E,x0,xf,u0,uf,p,t0,tf,data);
            if intres_reg
                if strcmp(data.resmin.dataNLP.options.derivatives,'adigator')
                    sol{phaseNo}.gradCost(1,ResNorm_vec.dY_location)=sol{phaseNo}.gradCost(1,ResNorm_vec.dY_location)+ResNorm_vec.dY';
                else
                    sol{phaseNo}.gradCost=sol{phaseNo}.gradCost+resMin_gradientCost_ModeMinRes(X,U,P,tau,t0,tf,data.resmin);
                end
            end
            solution=sol{phaseNo}.gradCost;
            

        case{'const'}
            if strcmp(data.options.derivatives,'adigator')
                sol{phaseNo}.const=[const_vec_Adigator.f;
                          avrc(X,U,P,t,data)';
                          b(x0,xf,u0,uf,p,t0,tf,vdat)];
            else
                X_vect=reshape(X',n*M,1);
                g_mat=g(X,U,P,t,vdat);
                g_vect=reshape(g_mat',M*ng,1);
                sol{phaseNo}.const=[(x0-data.x0t)*data.cx0;
                          mp.A*X_vect+mp.B*reshape((tf-t0)*f(X,U,P,t,vdat)',M*n,1);
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
            
        case{'costHessian'}
            
          if strcmp(data.options.derivatives,'analytic') && (~isfield(data.options,'forceNumericHes') || ~data.options.forceNumericHes)
                [Lzz,Ezz]=hessianAN(L,f,g,sol{phaseNo}.Jf,sol{phaseNo}.JL,X,U,P,tau,E,b,x0,...
                                                        xf,u0,uf,p,t0,tf,data);          
          elseif strcmp(data.options.derivatives,'adigator')
                [Lzz,Ezz]=hessianCDAdigator(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data); 
          else
                [Lzz,Ezz]=hessianCD(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data); 
          end
          solution=Lzz+Ezz;
          
        case{'hessian'}
            
          if strcmp(data.options.derivatives,'analytic') && (~isfield(data.options,'forceNumericHes') || ~data.options.forceNumericHes)
                [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,sol{phaseNo}.Jf,sol{phaseNo}.JL,X,U,P,tau,E,b,x0,...
                                                        xf,u0,uf,p,t0,tf,data);          
                fgzz=fzz+gzz;
          elseif strcmp(data.options.derivatives,'adigator')
              if intres_reg
                 data.ResNorm_vec=ResNorm_vec;
                 [Lzz,Ezz,fgzz,bzz,Reszz]=hessianCDAdigator(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data); 
              else
                [Lzz,Ezz,fgzz,bzz]=hessianCDAdigator(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data); 
              end
          else
              if intres_reg
                  [Lzz,Ezz,fzz,gzz,bzz,Reszz]=hessianCD(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data); 
              else
                  [Lzz,Ezz,fzz,gzz,bzz]=hessianCD(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data); 
              end
                fgzz=fzz+gzz;
          end
%           
          if intres_reg
              hessc=data.sigma*(Lzz+Ezz+Reszz)+fgzz+bzz;
          else
              hessc=data.sigma*(Lzz+Ezz)+fgzz+bzz;
          end
          hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
          sol{phaseNo}.hessian=tril(hessc);    
          
          if strcmp(data.options.NLPsolver,'builtinSQP')
            sol{phaseNo}.hessian=sol{phaseNo}.hessian+sol{phaseNo}.hessian'-diag([diag(sol{phaseNo}.hessian)]);
          end
          if strcmp(data.options.NLPsolver,'worhp')
              solution=sol{phaseNo}.hessian(sub2ind(size(sol{phaseNo}.hessian),data.jacStruct_HMRow,data.jacStruct_HMCol));
          else
              solution=sol{phaseNo}.hessian;
          end
       
       case{'constTest'}
            X_vect=reshape(X',n*M,1);
            g_mat=g(X,U,P,t,vdat);
            if strcmp(data.options.discretization,'discrete') || strcmp(data.options.discretization,'euler')
                g_mat(end-ng+1:end,:)=0;
            end
            g_vect=reshape(g_mat',M*ng,1);
            solution.fc=[(x0-data.x0t)*data.cx0;mp.A*X_vect+mp.B*reshape((tf-t0)*f(X,U,P,t,vdat)',M*n,1)];
            solution.gc=g_vect(data.gAllidx);
            solution.avrcc=avrc(X,U,P,t,data)';
            solution.bc=b(x0,xf,u0,uf,p,t0,tf,vdat);    
       case{'gradCostTest'}
            [solution.Lz,solution.Ez,~]=gradientCost(L,X,Xr,U,Ur,P,tau,E,x0,xf,u0,uf,p,t0,tf,data);
       case{'jacConstTest'}
            [solution.fzd,solution.gzd,solution.rcz,solution.bz]=jacobianFD(f,g,avrc,X,U,P,tau,b,x0,xf,u0,uf,p,t0,tf,data);
       case{'hessianTest'}
            [solution.Lzz,solution.Ezz]=hessianCostCD(L,f,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data); 

          
          
%           if strcmp(data.options.NLPsolver,'OSQP')
%               hessc = triu(hessc.',1) + tril(hessc);
%           end
            

    end

catch e
    fprintf(1,'There was an error! The message was:\n%s \n',e.message);
    error(['Error encountered when evaluating ' required ' function'])
end
