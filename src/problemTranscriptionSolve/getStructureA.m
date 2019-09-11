function  structure=getStructureA(data)

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
    else
      structure.dEdtf= 0;  
    end
    if (nt>=2)&&(~isempty(dE.dt0))
      structure.dEdt0 = 1;
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
else
    if (nt>=2)
        structure.dEdt0=spones(1*ones(1,1));
        structure.dEdtf=spones(1*ones(1,1));
    else
        structure.dEdtf=spones(nt*ones(1,nt));
    end
    structure.dEdp=spones(np*ones(1,np));
    structure.dEdu0=ones(1,m); 
    structure.dEdx0=ones(1,n);
    structure.dEduf=ones(1,m);
    structure.dEdxf=ones(1,n);
end    


% Preallocate memory
dLdx=ones(1,n);dfdx=ones(n);
dLdu=ones(1,m);dfdu=ones(n,m);

dLdt=ones(1,nt);

if ng
 dgdx=ones(ng,n);
 dgdu=ones(ng,m);
 if np
  dgdp=ones(ng,np);
 else
  dgdp=zeros(ng,np);  
 end
 dgdt=ones(ng,nt);
else
 dgdx=zeros(ng,n);
 dgdu=zeros(ng,m);
 dgdp=zeros(ng,np);
 dgdt=zeros(ng,nt);
end

if nrc
    if strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR')
        M=sum(npd);
        drcdx=ones(nrc/M,n);
        drcdu=ones(nrc/M,m);
        drcdp=ones(nrc/M,np);  
        drcdt=ones(nrc/M,nt);
    else
        drcdx=ones(nrc/((M-1)/ns),n);
        drcdu=ones(nrc/((M-1)/ns),m);
        drcdp=ones(nrc/((M-1)/ns),np);  
        drcdt=ones(nrc/((M-1)/ns),nt);
    end
else
    if strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR')
        M=sum(npd);
        drcdx=zeros(nrc/M,n);
        drcdu=zeros(nrc/M,m);
        drcdp=zeros(nrc/M,np);  
        drcdt=zeros(nrc/M,nt);
    else
        drcdx=zeros(nrc/((M-1)/ns),n);
        drcdu=zeros(nrc/((M-1)/ns),m);
        drcdp=zeros(nrc/((M-1)/ns),np);  
        drcdt=zeros(nrc/((M-1)/ns),nt);
    end
end



if nb
   %%%%%%%%%%%%%%%%%%%%%% 
    
     if db.flag==1
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
    else
        if (nt>=2)
           structure.dbdt0=spones(1*ones(nb,1));
           structure.dbdtf=spones(1*ones(nb,1));
        else
           structure.dbdtf=spones(nt*ones(nb,nt)); 
        end
        
        structure.dbdp=spones(np*ones(nb,np));
        structure.dbdu0=ones(nb,m); 
        structure.dbdx0=ones(nb,n);
        structure.dbduf=ones(nb,m);
        structure.dbdxf=ones(nb,n);
    end    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
      structure.dbdx0=zeros(nb,n);
      structure.dbdxf=zeros(nb,n);
      structure.dbdu0=zeros(nb,m);
      structure.dbduf=zeros(nb,m);
      structure.dbdp=zeros(nb,np);
      if (nt>=2)
        structure.dbdt0=zeros(nb,1);  
        structure.dbdtf=zeros(nb,1);
      else
        structure.dbdtf=zeros(nb,nt);
      end
        
end
    


if np
  dLdp=ones(1,np);
  dfdp=ones(n,np);
else  
  dLdp=zeros(1,np);
  dfdp=zeros(n,np);
end   
    


% 
%  Sturctures for the jacobians 
%---------------------------------------


structure.dfdx=spones(dfdx);
structure.dLdx=spones(dLdx);


structure.dfdu=spones(dfdu);
structure.dLdu=spones(dLdu);



  
structure.dfdp=spones(dfdp);structure.dLdp=spones(dLdp);
structure.dLdt=spones(dLdt);
structure.dgdt=spones(dgdt);

structure.drcdx=spones(drcdx);structure.drcdu=spones(drcdu);
structure.drcdp=spones(drcdp);structure.drcdt=spones(drcdt);



structure.dgdx=spones(dgdx);structure.dgdu=spones(dgdu);
structure.dgdp=spones(dgdp);structure.dgdt=spones(dgdt);



%------------- END OF CODE --------------

