function [varargout]=getPertubations(sparsity,sizes,data)
%getPertubations - Generate pertubation vectors 
%
% Syntax:  [vector,index]=getPertubations(sparsity,sizes,data)
%
% Inputs: Defined in transcribeOCP.m
%
% Outputs:
%    vector - structure containing the pertubation vectors
%    index  - structure containing the mapping to the jacobian
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
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,~]=deal(sizes{1:13});
nrc=nrcl+nrcu+nrce;
nz=nt+np+M*n+N*m;

% Get sparsity templates
dLdp=sparsity.dLdp;     dLdx=sparsity.dLdx;        dLdu=sparsity.dLdu;
dLdt=sparsity.dLdt;   

dfdp=sparsity.dfdp;        dfdx=sparsity.dfdx;
dfdu=sparsity.dfdu;     

drcdp=sparsity.drcdp;        drcdt=sparsity.drcdt;

dgdp=sparsity.dgdp;        dgdt=sparsity.dgdt;       
dbdx0=sparsity.dbdx0;   
dbdxf=sparsity.dbdxf;   dbdu0=sparsity.dbdu0;      dbduf=sparsity.dbduf;  
dbdp=sparsity.dbdp;     dbdt0=sparsity.dbdt0;      dbdtf=sparsity.dbdtf;



if  strcmp(data.options.derivatives,'numeric') || strcmp(data.options.derivatives,'adigator') 
  index.dE.flag=0;
  index.dL.flag=0;
  index.dR.flag=0;
  index.df.flag=0;
  index.dg.flag=0;
  index.drc.flag=0;
  index.db.flag=0;

else
 index.dE.flag=sparsity.dE.flag;
 index.dL.flag=sparsity.dL.flag;
 index.df.flag=sparsity.df.flag;
 index.dg.flag=sparsity.dg.flag;
 index.drc.flag=0;
 index.dR.flag=0;
 index.db.flag=sparsity.db.flag;
end


ez_scaled=[ones(nt,1);ones(np,1);ones(n,1);ones(m,1)];
ezf_scaled=[ones(nt,1);ones(np,1);ones(n,1);ones(m,1);ones(n,1);ones(m,1)];



% 1.Get pertubations and indices for dLdz. z:=[tf p x0 u0 .... xf uf]'

idx0=find([dLdt dLdp dLdx dLdu]);
nfd=nnz(idx0);
idx=zeros(M,nfd);

et0=cell(1,nfd);etf=cell(1,nfd);et=cell(1,nfd);ep=cell(1,nfd);
ex=cell(1,nfd);eu=cell(1,nfd);
e0=1:M;



for i=1:nfd                                    % Allocate memory 
    et0{i}=zeros(1);etf{i}=zeros(1);et{i}=spalloc(M,nt,M);ep{i}=spalloc(M,np,M);
    ex{i}=spalloc(M,n,M);eu{i}=spalloc(M,m,M);
end


for i=1:nfd

 if idx0(i)<=nt                                         % Pertubation of tf
%    etf{i}=1;idx(:,i)=1;
    et{i}=sparse(e0,idx0(i),1,M,nt);
    idx(:,i)=idx0(i);
    if nt==2 && idx0(i)==1
        et0{i}=1;
    elseif nt==2 && idx0(i)==2
        etf{i}=1;
    end
 end

 if (idx0(i)>nt)&&(idx0(i)<nt+np+1)                    % Pertubation of p
    ep{i}=sparse(e0,idx0(i)-nt,1,M,np);
    idx(:,i)=idx0(i);
 end

 
 if (idx0(i)>nt+np)&&(idx0(i)<nt+np+n+1)               % Pertubation of x
     ex{i}=sparse(e0,idx0(i)-nt-np,1,M,n);
     idx(:,i)=find(data.map.xV*repmat(sparse(idx0(i)-nt-np,1,1,n,1),M,1))';
 end


 if (idx0(i)>nt+np+n)&&(idx0(i)<nt+np+n+m+1)           % Pertubation of u
     eu{i}=sparse(e0,idx0(i)-nt-np-n,1,M,m);
     idx(:,i)=find(data.map.uV*repmat(sparse(idx0(i)-nt-np-n,1,1,m,1),N,1))';
 end
end

idx=sparse(idx);

vector.Ly.et0=et0;vector.Ly.etf=etf;vector.Ly.et=et;vector.Ly.ep=ep;
vector.Ly.ex=ex;vector.Ly.eu=eu;vector.Ly.ez=ez_scaled;
index.Ly=idx;
 



% 2. Get pertubations and indices for dEdy. y:=[tf p x0 u0  xf uf]'
dEdx0=sparsity.dEdx0; 
dEdxf=sparsity.dEdxf;
dEdu0=sparsity.dEdu0;
dEduf=sparsity.dEduf; 
dEdp=sparsity.dEdp;
dEdt0=sparsity.dEdt0;
dEdtf=sparsity.dEdtf;

 idx0=find([dEdt0 dEdtf dEdp dEdx0 dEdu0 dEdxf dEduf]);
 idx=idx0;nfd=nnz(idx);
 et0=zeros(1,nfd);etf=zeros(1,nfd);ep=zeros(np,nfd);                      % Allocate memory 
 ex0=zeros(n,nfd);eu0=zeros(m,nfd);
 exf=zeros(n,nfd);euf=zeros(m,nfd);

 for i=1:nfd
    
    if idx0(i)==nt                                     % Pertubation of tf
        if nt==1
            etf(i)=1;
        else
            if i~=1
                et0(i-1)=1;
            end
            etf(i)=1;
        end
    end
   if (idx0(i)>nt)&&(idx0(i)<nt+np+1)                 % Pertubation of p
        ep(:,i)=sparse(idx0(i)-nt,1,1,np,1);
    end

    if (idx0(i)>nt+np)&&(idx0(i)<nt+np+n+1)            % Pertubation of x0
        ex0(:,i)=sparse(idx0(i)-nt-np,1,1,n,1);
    end

    if (idx0(i)>nt+np+n)&&(idx0(i)<nt+np+n+m+1)        % Pertubation of u0
        eu0(:,i)=sparse(idx0(i)-nt-np-n,1,1,m,1);
        idx(i)=idx(i)+(M/N-1)*n;                     
    end
    if (idx0(i)>nt+np+n+m)&&(idx0(i)<nt+np+2*n+m+1)    % Pertubation of xf
        exf(:,i)=sparse(idx0(i)-nt-np-n-m,1,1,n,1);
        idx(i)=idx(i)+(M-2)*n+(N-2)*m;
    end
    if (idx0(i)>nt+np+2*n+m)&&(idx0(i)<nt+np+2*n+2*m+1)    % Pertubation of uf
        euf(:,i)=sparse(idx0(i)-nt-np-2*n-m,1,1,m,1);
        idx(i)=idx(i)+(M-2)*n+(N-2)*m;
    end
end
index.Ey=idx;
vector.Ey.et0=et0;vector.Ey.etf=etf;vector.Ey.ep=ep;
vector.Ey.ex0=ex0;vector.Ey.eu0=eu0;
vector.Ey.exf=exf;vector.Ey.euf=euf;vector.Ey.ez=ezf_scaled;
% For multiphase linkage constraint
if nargout==1 || nargout==3
    infoForLinkConst.x0idx=index.Ey(nt+np+1:nt+np+n)-nt-np;
    infoForLinkConst.u0idx=index.Ey(nt+np+n+1:nt+np+n+m)-nt-np;
    infoForLinkConst.xfidx=index.Ey(nt+np+n+m+1:nt+np+n*2+m)-nt-np;
    infoForLinkConst.ufidx=index.Ey(nt+np+n*2+m+1:nt+np+n*2+m*2)-nt-np;
    infoForLinkConst.nxu0f=n*2+m*2;
end

% dfdz where z=[tf p x0 ...u0... xf uf];
% -------------------------------------
Fxu=[kron(speye((M)/N),dfdx),repmat(dfdu,(M)/N,1)];
df=[repmat(ones(n,nt),M,1) repmat(dfdp,M,1) kron(speye(N),Fxu)];

colFxu=repmat(overlapping(Fxu),N,1);
if isempty(dfdp); dp=0;else dp=max(overlapping(dfdp));end


colgroup=[[1:nt]';ones(~(np==0))*(overlapping(dfdp)+nt);...      
         nt+dp+colFxu];ngcol=max(colgroup);                      


% eg=spalloc(nz,ngcol,nz);
eg=zeros(nz,ngcol);
vector.f.et0=zeros(1,ngcol);vector.f.et0(1)=1;
vector.f.etf=zeros(1,ngcol);vector.f.etf(2)=1;
vector.f.ep=eg;vector.f.ex=eg;vector.f.eu=eg;
ixf=zeros(n*M,ngcol);
% ixf=spalloc(n*M,ngcol,n*M*ngcol);


for i=1:ngcol
    eg(colgroup==i,i)=1;
    [a,b]=find(df*spdiags(eg(:,i),0,nz,nz));
    ixf(a,i)=b;
    ixf(ixf(:,i)==0,i)=1;
end


ex=data.map.Vx*eg;eu=data.map.Vu*eg;


vector.f.ex=cell(ngcol,1);vector.f.eu=cell(ngcol,1);vector.f.ep=cell(ngcol,1);vector.f.ez=ez_scaled;


for i=1:ngcol
    vector.f.ex{i}=reshape(ex(:,i)',n,M)';
    vector.f.eu{i}=kron(reshape(eu(:,i)',m,N)',ones(M/N,1)); 
    vector.f.ep{i}=repmat(eg(nt+1:nt+np,i)',M,1);
end

index.f=ixf;



% dgdz where z=[tf p x0 ...u0...  uf xf]; 

if ng
Gxu=[kron(speye(M/N),sparse(sparsity.dgdx)),repmat(sparse(sparsity.dgdu),M/N,1)];

dg=[repmat(dgdt,M,1) repmat(dgdp,M,1) kron(speye(N),Gxu)];

   
cGxu=overlapping(Gxu);   
colGxu=repmat(cGxu,N,1);

if isempty(dgdp); dp=0;else dp=max(overlapping(dgdp));end


colgroup=[[1:nt]';ones(~(np==0))*(overlapping(dgdp)+nt);...
         nt+dp+colGxu];ngcol=max(colgroup);


eg=zeros(nz,ngcol);
vector.g.et0=zeros(1,ngcol);vector.g.et0(1)=1;
vector.g.etf=zeros(1,ngcol);vector.g.etf(2)=1;
vector.g.ep=eg;vector.g.ex=eg;vector.g.eu=eg;
ig=zeros(ng*M,ngcol);


for i=1:ngcol
    eg(colgroup==i,i)=1;
    [a,b]=find(dg*spdiags(eg(:,i),0,nz,nz));
    
    ig(a,i)=b;
       
    ig(ig(:,i)==0,i)=1;
  
end
eg=sparse(eg);
ig=sparse(ig);

ex=data.map.Vx*eg;eu=data.map.Vu*eg;

vector.g.ex=cell(ngcol,1);vector.g.eu=cell(ngcol,1);vector.g.ep=cell(ngcol,1);vector.g.ez=ez_scaled;

for i=1:ngcol
    vector.g.ex{i}=reshape(ex(:,i)',n,M)';
    vector.g.eu{i}=kron(reshape(eu(:,i)',m,N)',ones((M)/N,1));
    vector.g.ep{i}=repmat(eg(nt+1:nt+np,i)',M,1);
end

index.g=ig;
end

% dgdz where z=[tf p x0 ...u0...  uf xf]; 


if nrc
RCxu=[kron(speye(M/N),sparse(sparsity.drcdx)),repmat(sparse(sparsity.drcdu),M/N,1)];

drc=[repmat(drcdt,M,1) repmat(drcdp,M,1) kron(speye(M),RCxu)];
   
cRCxu=overlapping(RCxu);   
colRCxu=repmat(cRCxu,N,1);
if isempty(drcdp); dp=0;else dp=max(overlapping(drcdp));end


colgroup=[[1:nt]';ones(~(np==0))*(overlapping(drcdp)+nt);...
         nt+dp+colRCxu];nrccol=max(colgroup);


erc=zeros(nz,nrccol);
vector.rc.et0=zeros(1,nrccol);vector.rc.et0(1)=1;
vector.rc.etf=zeros(1,nrccol);vector.rc.etf(2)=1;
vector.rc.ep=erc;vector.rc.ex=erc;vector.rc.eu=erc;
irc=zeros(nrc,nrccol);

for i=1:nrccol
    erc(colgroup==i,i)=1;
    [a,b]=find(drc*spdiags(erc(:,i),0,nz,nz));
    
    irc(a,i)=b;
       
    irc(irc(:,i)==0,i)=1;
  
end
erc=sparse(erc);
irc=sparse(irc);


ex=data.map.Vx*erc;eu=data.map.Vu*erc;

vector.rc.ex=cell(nrccol,1);vector.rc.eu=cell(nrccol,1);vector.rc.ep=cell(nrccol,1);vector.rc.ez=ez_scaled;

for i=1:nrccol
    vector.rc.ex{i}=reshape(ex(:,i)',n,M)';
    vector.rc.eu{i}=kron(reshape(eu(:,i)',m,N)',ones(M/N,1));
    vector.rc.ep{i}=repmat(erc(nt+1:nt+np,i)',M,1);
end

irc=irc(1:ns:end,:);
idxicrm=1:size(irc,1)/M;
irc(idxicrm*M,:)=[];
index.rc=irc;
end



% dbdzb where zb=[tf p x0 u0 uf xf];
% -----------------------------------



db=[dbdt0 dbdtf dbdp dbdx0 dbdu0 dbdxf dbduf];
colgroup=overlapping(db);nbcol=max(colgroup);



eb=zeros(2*(n+m)+np+nt,nbcol);
vector.b.et0=eb(1,:);
vector.b.etf=eb(1,:);
vector.b.ep=eb(1:np,:);
vector.b.ex0=eb(1:n,:);vector.b.eu0=eb(1:m,:);
vector.b.exf=eb(1:n,:);vector.b.euf=eb(1:m,:);

ib=zeros(nb,nbcol);

if nb
for i=1:nbcol
    eb(colgroup==i,i)=1;
    for j=1:nb
       ib(j,i)=intersect(find(eb(:,i)),find(db(j,:)));
    end
end

if nt
    vector.b.et0=eb(1,:);
    vector.b.etf=eb(2,:);
end
if np; vector.b.ep=eb(nt+1:nt+np,:); end

vector.b.ex0=eb(nt+np+1:nt+np+n,:);
vector.b.eu0=eb(nt+np+n+1:nt+np+n+m,:);
vector.b.exf=eb(nt+np+n+m+1:nt+np+2*n+m,:);
vector.b.euf=eb(nt+np+2*n+m+1:nt+np+2*n+2*m,:);
vector.b.ez=ezf_scaled;
end


% shift index
ibm=ib;
ibm((ib>nt+np+n)&(ib<nt+np+n+m+1))=ib((ib>nt+np+n)&(ib<nt+np+n+m+1))+((M)/N-1)*n;   
ibm(ib>nt+np+n+m)=nz-nt-np-2*n-2*m+ib((ib>nt+np+n+m));
index.idb=ib;
index.b=ibm;


%% Output
if nargout==1
    varargout{1}=infoForLinkConst;
elseif nargout==2
    disp('Generate finite-difference pertubation vectors')
    varargout{1}=vector;
    varargout{2}=index;
else
    disp('Generate finite-difference pertubation vectors')
    varargout{1}=vector;
    varargout{2}=index;
    varargout{3}=infoForLinkConst;
end



%------------- END OF CODE --------------
  