function  structure=getStructure(sizes,method)

%GETSTRUCTURE - Generate sparsity templates. 
%
% Syntax:  structure=getStucture(functions,sizes,method,data)
%
% Inputs: Defined in transcribeOCP.m
%
% Outputs:
%    structure - Structure containing the sparsity templates
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


disp('Determine sparsity structure');

% Get dimensions
if strcmp(method,'globalLGR') || strcmp(method,'hpLGR')
    [nt,np,n,m,ng,nb,M,~,ns,npd,~,~,~,nrcl,nrcu,nrce]=deal(sizes{:});
else
    [nt,np,n,m,ng,nb,M,~,ns,nrcl,nrcu,nrce]=deal(sizes{:});
end
nrc=nrcl+nrcu+nrce;

% Preallocate memory 
dLdx=ones(1,n);dEdx0=ones(1,n);dEdxf=ones(1,n);dfdx=ones(n);


dLdu=ones(1,m);dfdu=ones(n,m);dEdu0=ones(1,m);dEduf=ones(1,m);

if strcmp(method,'globalLGR') || strcmp(method,'hpLGR')
    dLdt=ones(1,nt);
else
    if (nt==1)                                      
        dLdtf=1;
    else
        dLdtf=[];
    end
end

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
    if strcmp(method,'globalLGR') || strcmp(method,'hpLGR')
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
    if strcmp(method,'globalLGR') || strcmp(method,'hpLGR')
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
	dbdx0=ones(nb,n);
	dbdxf=ones(nb,n);  
	dbdu0=ones(nb,m);
	dbduf=ones(nb,m);
	if np
        dbdp=ones(nb,np);  
    else
        dbdp=zeros(nb,np);  
    end
    if (nt==1)
        dbdtf=ones(nb,nt);
    elseif (nt>=2)
        dbdt0=ones(nb,1);
        dbdtf=ones(nb,1);
    else
        dbdtf=zeros(nb,nt);    
    end   
else
	dbdx0=zeros(nb,n);
	dbdxf=zeros(nb,n);
	dbdu0=zeros(nb,m);
	dbduf=zeros(nb,m);
	dbdp=zeros(nb,np);
	dbdt0=zeros(nb,nt);  
	dbdtf=zeros(nb,nt);  
end
    

if np
	dLdp=ones(1,np);
	dfdp=ones(n,np);
	dEdp=ones(1,np);
else  
	dLdp=zeros(1,np);
	dfdp=zeros(n,np);
	dEdp=zeros(1,np);
end   
    


% 
%  Sturctures for the jacobians 
%---------------------------------------

structure.dfdx=spones(dfdx);
structure.dLdx=spones(dLdx);

structure.dEdx0=spones(dEdx0);
structure.dEdxf=spones(dEdxf);
structure.dfdu=spones(dfdu);


structure.dLdu=spones(dLdu);

structure.dEdu0=spones(dEdu0);
structure.dEduf=spones(dEduf);

  
structure.dfdp=spones(dfdp);structure.dLdp=spones(dLdp);
structure.dEdp=spones(dEdp);


if strcmp(method,'globalLGR') || strcmp(method,'hpLGR')
    structure.dLdt=spones(dLdt);
else
    structure.dLdtf=spones(dLdtf);
end

if (nt==1)                                      % Get structure of E wrt tf
    dEdtf = 1;
elseif (nt>=2)                                      % Get structure of E wrt tf
    dEdt0 = 1;
    dEdtf = 1;
    structure.dEdt0=spones(dEdt0);
else
    dEdtf=[];
end

structure.dEdtf=spones(dEdtf);



structure.dgdx=spones(dgdx);structure.dgdu=spones(dgdu);
structure.dgdp=spones(dgdp);structure.dgdt=spones(dgdt);

structure.drcdx=spones(drcdx);structure.drcdu=spones(drcdu);
structure.drcdp=spones(drcdp);structure.drcdt=spones(drcdt);

structure.dbdx0=spones(dbdx0);structure.dbdxf=spones(dbdxf);
structure.dbdu0=spones(dbdu0);structure.dbduf=spones(dbduf);
structure.dbdp=spones(dbdp);structure.dbdtf=spones(dbdtf);

if nt>=2
    structure.dbdt0=spones(dbdt0);
end


%------------- END OF CODE --------------

