function  [structure, data]=getStructureA(structure, data)

%GETSTRUCTUREA - Generate sparsity templates when the analytic option has been selected
%
% Syntax:  structure=getStructureA(functions,sizes,data)
%
% Inputs: Defined in transcribeOCP.m
%
% Outputs:
%    structure - Structure containing the sparsity templates
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

% disp('Determine sparsity structure');

% Get dimensions
if strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR')
    [nt,np,n,m,ng,nb,M,~,ns,npd,~,~,~,nrcl,nrcu,nrce,~]=deal(data.sizes{1:17});
else
    [nt,np,n,m,ng,nb,M,~,ns,nrcl,nrcu,nrce,~]=deal(data.sizes{1:13});
end
nrc=nrcl+nrcu+nrce;


% Check the structure of the derivatives for the stage cost L and the terminal cost E

[L,E,f,g,~,b]=deal(data.functions{:});
gradCost=data.analyticDeriv.gradCost;
jacConst=data.analyticDeriv.jacConst;

[dL,dE]=gradCost(L,rand(1,n),rand(1,n),rand(1,m),rand(1,m),rand(1,np),0,E,rand(n,1),rand(n,1),rand(m,1),rand(m,1),rand(np,1),rand(1,1),rand(1,1),data);
[df,dg,db]=jacConst(f,g,rand(1,n),rand(1,m),rand(1,np),0,b,rand(n,1),rand(n,1),rand(m,1),rand(m,1),rand(np,1),rand(1,1),rand(1,1),data);
structure.dE.flag=dE.flag;
structure.dL.flag=dL.flag;
structure.df.flag=df.flag;
structure.dg.flag=dg.flag;
structure.db.flag=db.flag;



if dE.flag==1
    if (nt~=0)&&(~isempty(dE.dtf))
      structure.dEdtf = 1;
    elseif nt==0 && isempty(dE.dtf)
      structure.dEdtf=[];  
    else
      structure.dEdtf= 0;  
    end
    if (nt>=2)&&(~isempty(dE.dt0))
      structure.dEdt0 = 1;
    elseif nt==0 && isempty(dE.dtf)
      structure.dEdt0=[];  
    else
      structure.dEdt0= 0;  
    end
    if np&&(~isempty(dE.dp))
      structure.dEdp=ones(1,np);
    else   
      structure.dEdp=zeros(1,np);  
    end
    if ~isempty(dE.dx0)
      structure.dEdx0=ones(1,n);
    else
      structure.dEdx0=zeros(1,n);  
    end
    if ~isempty(dE.du0)
      structure.dEdu0=ones(1,m);
    else
      structure.dEdu0=zeros(1,m);
    end
    if ~isempty(dE.duf)
     structure.dEduf=ones(1,m);
    else 
      structure.dEduf=zeros(1,m);
    end    
    if ~isempty(dE.dxf)  
     structure.dEdxf=ones(1,n);
    else
     structure.dEdxf=zeros(1,n);
    end
end    

data.FD.vindex.E=[];
% if dE.flag==1 
    if strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR')
        if dE.flag==0 || ~isempty(dE.dx0)
          data.FD.vindex.E = [data.FD.vindex.E, 1:n];
        end
        if dE.flag==0 || ~isempty(dE.dxf)  
          data.FD.vindex.E = [data.FD.vindex.E, n+1:2*n];
        end
        if dE.flag==0 || ~isempty(dE.du0)
          data.FD.vindex.E = [data.FD.vindex.E, n*2+1:2*n+m];
        end
        if dE.flag==0 || ~isempty(dE.duf)
          data.FD.vindex.E = [data.FD.vindex.E, n*2+m+1:2*n+2*m];
        end    
        if np&&(dE.flag==0 || ~isempty(dE.dp))
          data.FD.vindex.E = [data.FD.vindex.E, n*2+m*2+1:2*n+2*m+np];
        end
        if (nt>=2)&&(dE.flag==0 || ~isempty(dE.dt0))
          data.FD.vindex.E = [data.FD.vindex.E, n*2+m*2+np+1];
        end
        if (nt~=0)&&(dE.flag==0 || ~isempty(dE.dtf))
          data.FD.vindex.E = [data.FD.vindex.E, n*2+m*2+np+1+1];
        end
    else
        if (nt~=0)&&(dE.flag==0 || ~isempty(dE.dtf))
          data.FD.vindex.E = [data.FD.vindex.E, 1];
        end
        if (nt>=2)&&(dE.flag==0 || ~isempty(dE.dt0))
          data.FD.vindex.E = [data.FD.vindex.E, 2];
        end
        if np&&(dE.flag==0 || ~isempty(dE.dp))
          data.FD.vindex.E = [data.FD.vindex.E, nt+1:nt+np];
        end
        if dE.flag==0 || ~isempty(dE.dx0)
          data.FD.vindex.E = [data.FD.vindex.E, nt+np+1:nt+np+n];
        end
        if dE.flag==0 || ~isempty(dE.du0)
          data.FD.vindex.E = [data.FD.vindex.E, nt+np+n+1:nt+np+n+m];
        end
        if dE.flag==0 || ~isempty(dE.dxf)  
          data.FD.vindex.E = [data.FD.vindex.E, nt+np+n+m+1:nt+np+n*2+m];
        end
        if dE.flag==0 || ~isempty(dE.duf)
          data.FD.vindex.E = [data.FD.vindex.E, nt+np+n*2+m+1:nt+np+n*2+m*2];
        end   
    end
% end    

% Preallocate memory
if dL.flag
    dLdx=ones(1,n);
    dLdu=ones(1,m);
    dLdt=ones(1,nt);
    if np
      dLdp=ones(1,np);
    else  
      dLdp=zeros(1,np);
    end   
    structure.dLdx=spones(dLdx);
    structure.dLdu=spones(dLdu);
    structure.dLdp=spones(dLdp);
    structure.dLdt=spones(dLdt);
end

if df.flag
    dfdx=ones(n);
    dfdu=ones(n,m);
    if np
      dfdp=ones(n,np);
    else  
      dfdp=zeros(n,np);
    end   
    structure.dfdx=spones(dfdx);
    structure.dfdu=spones(dfdu);
    structure.dfdp=spones(dfdp);
end

if ng && dg.flag
 dgdx=ones(ng,n);
 dgdu=ones(ng,m);
 if np
  dgdp=ones(ng,np);
 else
  dgdp=zeros(ng,np);  
 end
 dgdt=ones(ng,nt);
 structure.dgdx=spones(dgdx);structure.dgdu=spones(dgdu);
 structure.dgdp=spones(dgdp);structure.dgdt=spones(dgdt);
 structure.dgdt=spones(dgdt);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nb && db.flag==1
      if (nt==1)&&(~isempty(db.dtf))
        structure.dbdtf =ones(nb,1);
      elseif (nt>=2)
          if (~isempty(db.dtf))&&(~isempty(db.dt0))
                structure.dbdt0 =ones(nb,1);
                structure.dbdtf =ones(nb,1);
          elseif (isempty(db.dtf))&&(~isempty(db.dt0))
                structure.dbdt0 =ones(nb,1);
                structure.dbdtf =zeros(nb,1);
          elseif (~isempty(db.dtf))&&(isempty(db.dt0))
                structure.dbdt0 =zeros(nb,1);
                structure.dbdtf =ones(nb,1);
          else
                structure.dbdt0 =zeros(nb,1);
                structure.dbdtf =zeros(nb,1);
          end
      else
        structure.dbdtf=zeros(nb,nt);  
      end
      if np&&(~isempty(db.dp))
         structure.dbdp=ones(nb,np);
      else   
         structure.dbdp=zeros(nb,np);  
      end
      if ~isempty(db.dx0)
        structure.dbdx0=ones(nb,n);
      else
        structure.dbdx0=zeros(nb,n);  
      end
      if ~isempty(db.du0)
        structure.dbdu0=ones(nb,m);
      else
        structure.dbdu0=zeros(nb,m);
      end
      if ~isempty(db.duf)
       structure.dbduf=ones(nb,m);
      else 
        structure.dbduf=zeros(nb,m);
      end    
      if ~isempty(db.dxf)  
       structure.dbdxf=ones(nb,n);
      else
       structure.dbdxf=zeros(nb,n);
      end
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%



    


% 
%  Sturctures for the jacobians 
%---------------------------------------
















%------------- END OF CODE --------------

