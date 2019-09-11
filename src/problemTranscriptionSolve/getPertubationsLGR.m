function [varargout]=getPertubationsLGR(sparsity,sizes,data)

%getPertubationsLGR - Generate pertubation vectors for LGR collocation
%
% Syntax:  [vector,index,Jac_templete]=getPertubationsLGR(sparsity,sizes,data)
%
% Inputs: Defined in transcribeOCP.m
%
% Outputs:
%    vector - structure containing the pertubation vectors
%    index  - structure containing the mapping to the jacobian
%    Jac_templete - a Jacobian matrix template
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

% Get dimensions
[nt,np,n,m,ng,nb,M,N,~,~,~,~,~,nrcl,nrcu,nrce,~]=deal(sizes{1:17});
nrc=nrcl+nrcu+nrce;
nz=nt+np+(M+1)*n+M*m;

% Get sparsity templates
dLdp=sparsity.dLdp;     dLdx=sparsity.dLdx;        dLdu=sparsity.dLdu;
dLdt=sparsity.dLdt;   
                        
dfdp=sparsity.dfdp;        dfdx=sparsity.dfdx;
dfdu=sparsity.dfdu;     

dgdp=sparsity.dgdp;     dgdt=sparsity.dgdt;
dgdx=sparsity.dgdx;     dgdu=sparsity.dgdu;        

drcdp=sparsity.drcdp;     drcdt=sparsity.drcdt;
drcdx=sparsity.drcdx;     drcdu=sparsity.drcdu; 

dbdx0=sparsity.dbdx0;   
dbdxf=sparsity.dbdxf;   dbdu0=sparsity.dbdu0;      dbduf=sparsity.dbduf;  
dbdp=sparsity.dbdp;     dbdtf=sparsity.dbdtf;      dbdt0=sparsity.dbdt0;

if  strcmp(data.options.derivatives,'analytic') 
     index.dE.flag=sparsity.dE.flag;
     index.dL.flag=sparsity.dL.flag;
     index.df.flag=sparsity.df.flag;
     index.dg.flag=sparsity.dg.flag;
     index.db.flag=sparsity.db.flag;
else
      index.dE.flag=0;
      index.dL.flag=0;
      index.df.flag=0;
      index.dg.flag=0;
      index.db.flag=0;
end


% 1.Get pertubations and indices for dLdz. z:=[x0 .... xf u0 ... uf p t0 tf]' 

idx0{1}=find(dLdx);
idx0{2}=find(dLdu);
idx0{3}=find(dLdp);
idx0{4}=find(dLdt);
nfd=nnz(cell2mat(idx0));
idx=zeros(M,nfd);

ex=cell(1,nfd);eu=cell(1,nfd);
ep=cell(1,nfd);et=cell(1,nfd);et0=cell(1,nfd);etf=cell(1,nfd);
e0=1:M;



for i=1:nfd                                 % Allocate memory 
    
    ex{i}=spalloc(M,n,M);eu{i}=spalloc(M,m,M);
    ep{i}=spalloc(M,np,M);et{i}=zeros(1,nt);et0{i}=zeros(1);etf{i}=zeros(1);
end

i=1;
for k=1:size(idx0,2)
    
 if k==1               % Pertubation of x
     for local_i=1:size(idx0{1, k},2)
         ex{i}=sparse(e0,local_i,1,M,n);
         idx(:,i)=((local_i-1)*M+local_i):((M+1)*local_i-1);
         i=i+1;
     end
 end
 
 if k==2           % Pertubation of u
     for local_i=1:size(idx0{1, k},2)
         eu{i}=sparse(e0,local_i,1,M,m);
         idx(:,i)=((M+1)*n+1+(local_i-1)*M):((M+1)*n+local_i*M);
         i=i+1;
     end
 end
 
%  idx
 if k==3                    % Pertubation of p
    for local_i=1:size(idx0{1, k},2)
        ep{i}=sparse(e0,local_i,1,M,np);
        idx(:,i)=((M+1)*n+M*m+1+(local_i-1)):((M+1)*n+M*m+local_i);
        i=i+1;
    end
 end
%  idx
 if k==4                                         % Pertubation of t0
     

     for local_i=1:size(idx0{1, k},2)
         et{i}(local_i)=1;
         if local_i==1
             et0{i}=1;
         elseif local_i==size(idx0{1, k},2)
             etf{i}=1;
         end
         idx(:,i)=nz-nt+local_i;
         i=i+1;
     end
     

 end


 
end
idx=sparse(idx);
vector.Ly.et=et;
vector.Ly.et0=et0;vector.Ly.etf=etf;
vector.Ly.ep=ep;
vector.Ly.ex=ex;vector.Ly.eu=eu;
index.Ly=idx;



% 2. Get pertubations and indices for dEdy. y:=[x0 u0 uf xf p t0 tf ]'


    
dEdx0=sparsity.dEdx0; 
dEdxf=sparsity.dEdxf;
dEdu0=sparsity.dEdu0;
dEduf=sparsity.dEduf; 
dEdp=sparsity.dEdp;
dEdt0=sparsity.dEdt0;
dEdtf=sparsity.dEdtf;

idx0{1}=find(dEdx0);
idx0{2}=find(dEdxf);
idx0{3}=find(dEdu0);
idx0{4}=find(dEduf);
idx0{5}=find(dEdp);
idx0{6}=find(dEdt0);
idx0{7}=find(dEdtf);
nfd=nnz(cell2mat(idx0));

%  idx0=find([  dEdx0 dEdxf dEdu0  dEduf  dEdp dEdt0 dEdtf]);
 idx=zeros(1,nfd);
 et0=zeros(1,nfd);etf=zeros(1,nfd);ep=zeros(np,nfd);                      % Allocate memory 
 ex0=zeros(n,nfd);eu0=zeros(m,nfd);
 exf=zeros(n,nfd);euf=zeros(m,nfd);
 
i=1;
for k=1:size(idx0,2)

    if k==1               % Pertubation of x0
         for local_i=1:size(idx0{1, k},2)
            ex0(:,i)=sparse(local_i,1,1,n,1);
            idx(i)=(local_i-1)*(M+1)+1;
            i=i+1;
         end
    end
    if k==2               % Pertubation of xf
         for local_i=1:size(idx0{1, k},2)
            exf(:,i)=sparse(local_i,1,1,n,1);
            idx(i)=(local_i-1)*(M+1)+(M+1);
            i=i+1;
         end
    end
    if k==3               % Pertubation of u0
         for local_i=1:size(idx0{1, k},2)
            eu0(:,i)=sparse(local_i,1,1,m,1);
            idx(i)=n*(M+1)+(local_i-1)*M+1;
            i=i+1;
         end
    end
    if k==4               % Pertubation of uf
         for local_i=1:size(idx0{1, k},2)
            euf(:,i)=sparse(local_i,1,1,m,1);
            idx(i)=n*(M+1)+(local_i-1)*M+M;
            i=i+1;
         end
    end
    if k==5               % Pertubation of p
         for local_i=1:size(idx0{1, k},2)
            ep(:,i)=sparse(local_i,1,1,np,1);
            idx(i)=n*(M+1)+M*m+local_i;
            i=i+1;
         end
    end
    if k==6               % Pertubation of t0
         if ~isempty(idx0{k})
             et0(i)=1;
             idx(i)=nz-nt+1;
             i=i+1;
         end
    end
    if k==7               % Pertubation of t0
         if ~isempty(idx0{k})
             etf(i)=1;
             idx(i)=nz;
             i=i+1;
         end
    end
end
index.Ey=idx;
vector.Ey.et0=et0;vector.Ey.etf=etf;vector.Ey.ep=ep;
vector.Ey.ex0=ex0;vector.Ey.eu0=eu0;
vector.Ey.exf=exf;vector.Ey.euf=euf;
% For multiphase linkage constraint
if nargout==1 || nargout==4
    infoForLinkConst.x0idx=index.Ey(1:n);
    infoForLinkConst.u0idx=index.Ey(n*2+1:n*2+m);
    infoForLinkConst.xfidx=index.Ey(n+1:n*2);
    infoForLinkConst.ufidx=index.Ey(n*2+m+1:n*2+m*2);
    infoForLinkConst.nxu0f=n*2+m*2;
end

% dfdz where z=[ x0 ...xf u0...  uf p t0 tf];
% -------------------------------------

dfxu=sparse(M*n,(M+1)*n+M*m);
for fi=1:size(dfdx,1)
    for xi=1:size(dfdx,2)
        dfxu(((fi-1)*M+1):fi*M,((xi-1)*(M+1)+1):(xi*(M+1)-1))=kron(speye(M),dfdx(fi,xi));
    end
    for ui=1:size(dfdu,2)
        dfxu(((fi-1)*M+1):fi*M,(n*(M+1)+(ui-1)*M+1):(n*(M+1)+ui*M))=kron(speye(M),dfdu(fi,ui));
    end
end

colFxu=overlapping(dfxu);
if isempty(dfdp)
    dp=0;
else
    dp=max(overlapping(dfdp));
end


colgroup=[colFxu;ones(~(np==0))*(dp+max(colFxu));((1:nt)+max(colFxu)+dp)'];
ngcol=max(colgroup);    
     
eg=zeros(nz,ngcol);
vector.f.et0=zeros(1,ngcol);vector.f.et0(ngcol-nt+1)=1;
vector.f.etf=zeros(1,ngcol);vector.f.etf(ngcol)=1;
vector.f.ep=eg;vector.f.ex=eg;vector.f.eu=eg;vector.f.et=eg;
ixf=zeros(n*M,ngcol);


for i=1:ngcol
    eg(colgroup==i,i)=1;
    [a,b]=find([dfxu sparse(M*n,np+nt)]*spdiags(eg(:,i),0,nz,nz));
    ixf(a,i)=b;
    
    if i~=1
        ixf(ixf(:,i)==0,i)=ixf(end,i-1)+1;
    else
        ixf(ixf(:,i)==0,i)=1;
    end
end

for i=1:nt
    ixf(:,end-nt+i)=nz-nt+i;
end
ixf=sparse(ixf);
eg=sparse(eg);

ex=data.map.Vx*eg;eu=data.map.Vu*eg;et=data.map.Vt*eg;




vector.f.ex=cell(ngcol,1);vector.f.eu=cell(ngcol,1);vector.f.ep=cell(ngcol,1);vector.f.et=cell(nt,1);


for i=1:ngcol
   
    extemp=reshape(ex(:,i)',(M+1),n);
    vector.f.ex{i}=extemp(1:end-1,:);
    eus=kron(reshape(eu(:,i)',M,m),ones((M-1)/N,1)); 
    vector.f.eu{i}=eus;
    vector.f.ep{i}=repmat(eg((M+1)*n+M*m+1:(M+1)*n+M*m+np,i)',M,1);
    if i>ngcol-nt
        vector.f.et{i}=et(i-ngcol+nt,end-nt+1:end);
    else
        vector.f.et{i}=zeros(1,nt);
    end
end

index.f=ixf;
    
% dgdz where z=[ x0 ...xf u0...  uf p t0 ... tf]; 

if ng
    
    dgxu=sparse(M*ng,(M+1)*n+M*m);
    for gi=1:size(dgdx,1)
        for xi=1:size(dgdx,2)
            dgxu(((gi-1)*M+1):gi*M,((xi-1)*(M+1)+1):(xi*(M+1)-1))=kron(speye(M),dgdx(gi,xi));
        end
        for ui=1:size(dgdu,2)
            dgxu(((gi-1)*M+1):gi*M,(n*(M+1)+(ui-1)*M+1):(n*(M+1)+ui*M))=kron(speye(M),dgdu(gi,ui));
        end
    end


    colGxu=overlapping(dgxu);   
    
    if isempty(dgdp) 
        dp=0;
    else
        dp=max(overlapping(dgdp));
    end

    colgroup=[colGxu;ones(~(np==0))*(dp+max(colGxu));((1:nt)+max(colGxu)+dp)'];
    ngcol=max(colgroup);

    eg=zeros(nz,ngcol);
    vector.g.et0=zeros(1,ngcol);vector.g.et0(ngcol-nt+1)=1;
    vector.g.etf=zeros(1,ngcol);vector.g.etf(ngcol)=1;
    vector.g.ep=eg;vector.g.ex=eg;vector.g.eu=eg;vector.g.et=eg;
    ig=zeros(ng*M,ngcol);

    for i=1:ngcol
        eg(colgroup==i,i)=1;
        [a,b]=find([dgxu sparse(M*ng,np+nt)]*spdiags(eg(:,i),0,nz,nz));

        ig(a,i)=b;
        
        if i~=1
            ig(ig(:,i)==0,i)=ig(end,i-1)+1;
        else
            ig(ig(:,i)==0,i)=1;
        end


    end
    for i=1:nt
        ig(:,end-nt+i)=nz-nt+i;
    end
    eg=sparse(eg);
    ig=sparse(ig);
    ex=data.map.Vx*eg;eu=data.map.Vu*eg;et=data.map.Vt*eg;

    vector.g.ex=cell(ngcol,1);vector.g.eu=cell(ngcol,1);vector.g.ep=cell(ngcol,1);vector.g.et=cell(nt,1);

    
    for i=1:ngcol
        extemp=reshape(ex(:,i)',M+1,n);
        vector.g.ex{i}=extemp(1:end-1,:);
        ugs=kron(reshape(eu(:,i)',M,m),ones((M-1)/N,1));
        vector.g.eu{i}=ugs; % can be optimized
        vector.g.ep{i}=repmat(eg(end-nt-np+1:end-nt,i)',M,1);
        if i>ngcol-nt
            vector.g.et{i}=et(i-ngcol+nt,end-nt+1:end);
        else
            vector.g.et{i}=zeros(1,nt);
        end
    end

    index.g=ig;
end

%% drcdz where z=[ x0 ...xf u0...  uf p t0 ... tf]; 

if nrc
    
    drcxu=zeros((M+1)*nrc/M,(M+1)*n+M*m);
    for rci=1:size(drcdx,1)
        for xi=1:size(drcdx,2)
            drcxu(((rci-1)*(M+1)+1):rci*(M+1),((xi-1)*(M+1)+1):(xi*(M+1)))=kron(speye(M+1),drcdx(rci,xi));
        end
        for ui=1:size(dgdu,2)
            D_block_U=kron(speye(M+1),drcdu(rci,ui));
            drcxu(((rci-1)*(M+1)+1):rci*(M+1),(n*(M+1)+(ui-1)*(M)+1):(n*(M+1)+ui*(M)))=D_block_U(:,1:end-1);
        end
    end
    drcxu=sparse(drcxu);

    colRCxu=overlapping(drcxu);   
    
    if isempty(drcdp) 
        dp=0;
    else
        dp=max(overlapping(drcdp));
    end

    colgroup=[colRCxu;ones(~(np==0))*(dp+max(colRCxu));((1:nt)+max(colRCxu)+dp)'];
    nrccol=max(colgroup);

    erc=zeros(nz,nrccol);
    vector.rc.et0=zeros(1,nrccol);vector.rc.et0(nrccol-nt+1)=1;
    vector.rc.etf=zeros(1,nrccol);vector.rc.etf(nrccol)=1;
    vector.rc.ep=erc;vector.rc.ex=erc;vector.rc.eu=erc;vector.rc.et=erc;

    irc=zeros(nrc/M*(M+1),nrccol);

    for i=1:nrccol
        erc(colgroup==i,i)=1;
        [a,b]=find([drcxu sparse((M+1)*nrc/M,np+nt)]*spdiags(erc(:,i),0,nz,nz));

        irc(a,i)=b;
        
        if i~=1
            irc(irc(:,i)==0,i)=max(irc(1:end,i-1))+1;
        else
            irc(irc(:,i)==0,i)=1;
        end
    end

    for i=1:nt
        irc(:,end-nt+i)=nz-nt+i;
    end
    irc=sparse(irc);
    erc=sparse(erc);

    ex=data.map.Vx*erc;eu=data.map.Vu*erc;et=data.map.Vt*erc;

    vector.rc.ex=cell(nrccol,1);vector.rc.eu=cell(nrccol,1);vector.rc.ep=cell(nrccol,1);vector.rc.et=cell(nt,1);

    
    for i=1:ngcol
        extemp=reshape(ex(:,i)',M+1,n);
        vector.rc.ex{i}=extemp;
        ugs=reshape(eu(:,i)',M,m);
        vector.rc.eu{i}=ugs; % can be optimized
        vector.rc.ep{i}=repmat(eg(end-nt-np+1:end-nt,i)',M,1);
        if i>nrccol-nt
            vector.rc.et{i}=et(i-nrccol+nt,end-nt+1:end);
        else
            vector.rc.et{i}=zeros(1,nt);
        end
    end

    idxicrm=1:size(irc,1)/(M+1);
    irc(idxicrm*(M+1),:)=[];
    index.rc=irc;
end

% dbdzb where zb=[x0 xf u0 uf p t0 tf];
% -----------------------------------

    idxx0=1:M+1:(M+1)*(n-1)+1;
    idxxf=(M+1):(M+1):(M+1)*n;
    idxu0=(M+1)*n+1:(M):((M+1)*n+(M)*(m-1)+1);
    idxuf=(M+1)*n+(M):(M):((M+1)*n+(M)*m);
    idxp=zeros(1,np);
    for i=1:np
        idxp(i)=(M+1)*n+M*m+i;
    end
    idxt0=nz-nt+1;
    idxtf=nz;
    idxt=nz-nt+1:nz;
    
    if data.options.adaptseg==1 
        db=[dbdx0 dbdxf dbdu0 dbduf dbdp spones(ones(nb,nt))];
        ibms=[idxx0 idxxf idxu0 idxuf idxp idxt];
    else
        db=[dbdx0 dbdxf dbdu0 dbduf dbdp dbdt0 dbdtf];
        ibms=[idxx0 idxxf idxu0 idxuf idxp idxt0 idxtf];
    end

    colgroup=overlapping(db);nbcol=max(colgroup);

    eb=zeros(2*(n+m)+np+nt,nbcol);
    vector.b.et0=eb(1,:);
    vector.b.etf=eb(1,:);vector.b.ep=eb(1:np,:);vector.b.et=eb(1:nt,:);
    vector.b.ex0=eb(1:n,:);vector.b.eu0=eb(1:m,:);
    vector.b.exf=eb(1:n,:);vector.b.euf=eb(1:m,:);

    ib=zeros(nb,nbcol);


    for i=1:nbcol
        eb(colgroup==i,i)=1;
        for j=1:nb
           ib(j,i)=intersect(find(eb(:,i)),find(db(j,:)));
        end
    end



    vector.b.ex0=eb(1:n,:);
    vector.b.exf=eb(n+1:2*n,:);
    vector.b.eu0=eb(2*n+1:2*n+m,:);
    vector.b.euf=eb(2*n+m+1:2*n+2*m,:);
    vector.b.et=zeros(nbcol,nt);

    if np; vector.b.ep=eb(2*n+2*m+1:2*n+2*m+np,:); end
    if nt 
        if data.options.adaptseg==1 
            vector.b.et(nbcol-nt+1:nbcol,1:nt)=eb(end-nt+1:end,end-nt+1:end);
            vector.b.et0=eb(2*n+2*m+np+1,:);
            vector.b.etf=eb(2*n+2*m+np+nt,:); 
        else
            vector.b.et0=eb(2*n+2*m+np+1,:);
            vector.b.etf=eb(2*n+2*m+np+2,:); 
        end
    end

%%
if nb
    ibm=repmat(ibms(logical(db(1,:))),nb,1);
else
    ibm=repmat(ibms,nb,1);
end
index.idb=ib;
index.b=ibm;   


%% Save also the Jacobian Templete
Jac_templete.dfxu=dfxu;

if ng
    Jac_templete.dgxu=dgxu;
end


%% Output
if nargout==1
    varargout{1}=infoForLinkConst;
elseif nargout==3
    disp('Generate finite-difference pertubation vectors')
    varargout{1}=vector;
    varargout{2}=index;
    varargout{3}=Jac_templete;
else
    disp('Generate finite-difference pertubation vectors')
    varargout{1}=vector;
    varargout{2}=index;
    varargout{3}=Jac_templete;
    varargout{4}=infoForLinkConst;
end

%------------- END OF CODE --------------
  