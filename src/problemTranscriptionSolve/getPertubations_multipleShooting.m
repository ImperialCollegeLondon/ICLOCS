function [vector,index]=getPertubations_multipleShooting(sparsity,sizes,data)
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
disp('Generate finite-difference pertubation vectors')
% Get dimensions
[nt,np,n,m,ng,nb,M,N]=deal(sizes{1:8});

nz=nt+np+M*n+N*m;


% Get sparsity templates
dLdp=sparsity.dLdp;     dLdx=sparsity.dLdx;        dLdu=sparsity.dLdu;
dLdtf=sparsity.dLdtf;   
        
                        
dfdp=sparsity.dfdp;        dfdx=sparsity.dfdx;
dfdu=sparsity.dfdu;     dgdp=sparsity.dgdp;        dgdt=sparsity.dgdt;
dgdx=sparsity.dgdx;     dgdu=sparsity.dgdu;        dbdx0=sparsity.dbdx0;   
dbdxf=sparsity.dbdxf;   dbdu0=sparsity.dbdu0;      dbduf=sparsity.dbduf;  
dbdp=sparsity.dbdp;     dbdtf=sparsity.dbdtf;



if  strcmp(data.options.derivatives,'numeric') 
  index.dE.flag=0;
  index.dL.flag=0;
  index.df.flag=0;
  index.dg.flag=0;
  index.db.flag=0;

else
 index.dE.flag=sparsity.dE.flag;
 index.dL.flag=sparsity.dL.flag;
 index.df.flag=sparsity.df.flag;
 index.dg.flag=sparsity.dg.flag;
 index.db.flag=sparsity.db.flag;
end


% 1.Get pertubations and indices for dLdz. z:=[tf p x0 u0 .... xf uf]'

idx0=find([dLdtf dLdp dLdx dLdu]);
nfd=nnz(idx0);
idx=spalloc(M,nfd,M*nfd);

etf=cell(1,nfd);ep=cell(1,nfd);
ex=cell(1,nfd);eu=cell(1,nfd);
e0=1:M;



for i=1:nfd;                                    % Allocate memory 
    etf{i}=zeros(1);ep{i}=spalloc(M,np,M);
    ex{i}=spalloc(M,n,M);eu{i}=spalloc(M,m,M);
end


for i=1:nfd

 if idx0(i)==nt                                         % Pertubation of tf
   etf{i}=1;idx(:,i)=1;
 end

 if (idx0(i)>nt)&&(idx0(i)<nt+np+1);                    % Pertubation of p
    ep{i}=sparse(e0,idx0(i)-nt,1,M,np);
    idx(:,i)=idx0(i);
 end

 
 if (idx0(i)>nt+np)&&(idx0(i)<nt+np+n+1);               % Pertubation of x
     ex{i}=sparse(e0,idx0(i)-nt-np,1,M,n);
     idx(:,i)=find(data.map.xV*repmat(sparse(idx0(i)-nt-np,1,1,n,1),M,1))';
 end


 if (idx0(i)>nt+np+n)&&(idx0(i)<nt+np+n+m+1);           % Pertubation of u
     eu{i}=sparse(e0,idx0(i)-nt-np-n,1,M,m);
     idxs=kron(find(data.map.uV*repmat(sparse(idx0(i)-nt-np-n,1,1,m,1),N,1))',ones(1,(M-1)/N));
     idx(:,i)=[idxs, idxs(end)] ;
 end
 
 
end

vector.Ly.etf=etf;vector.Ly.ep=ep;
vector.Ly.ex=ex;vector.Ly.eu=eu;
index.Ly=idx;
 



% 2. Get pertubations and indices for dEdy. y:=[tf p x0 u0 uf xf]'


    
dEdx0=sparsity.dEdx0; 
dEdxf=sparsity.dEdxf;
dEdu0=sparsity.dEdu0;
dEduf=sparsity.dEduf; 
dEdp=sparsity.dEdp;
dEdtf=sparsity.dEdtf;

 idx0=find([dEdtf dEdp dEdx0 dEdu0  dEduf dEdxf]);
 idx=idx0;nfd=nnz(idx);
 etf=zeros(1,nfd);ep=zeros(np,nfd);                      % Allocate memory 
 ex0=zeros(n,nfd);eu0=zeros(m,nfd);
 exf=zeros(n,nfd);euf=zeros(m,nfd);

 for i=1:nfd
    
    if idx0(i)==nt;                                     % Pertubation of tf
        etf(i)=1;
    end
   if (idx0(i)>nt)&&(idx0(i)<nt+np+1);                 % Pertubation of p
        ep(:,i)=sparse(idx0(i)-nt,1,1,np,1);
    end

    if (idx0(i)>nt+np)&&(idx0(i)<nt+np+n+1);            % Pertubation of x0
        ex0(:,i)=sparse(idx0(i)-nt-np,1,1,n,1);
    end

    if (idx0(i)>nt+np+n)&&(idx0(i)<nt+np+n+m+1);        % Pertubation of u0
        eu0(:,i)=sparse(idx0(i)-nt-np-n,1,1,m,1);
        idx(i)=idx(i)+((M-1)/N-1)*n;                     
    end
    if (idx0(i)>nt+np+n+m)&&(idx0(i)<nt+np+n+2*m+1);    % Pertubation of uf
        euf(:,i)=sparse(idx0(i)-nt-np-n-m,1,1,m,1);
        idx(i)=idx(i)+(M-2)*n+(N-2)*m;
    end
    if (idx0(i)>nt+np+n+2*m)&&(idx0(i)<nt+np+2*n+2*m+1);    % Pertubation of xf
        exf(:,i)=sparse(idx0(i)-nt-np-n-2*m,1,1,n,1);
        idx(i)=idx(i)+(M-2)*n+(N-2)*m;
    end
end
index.Ey=idx;
vector.Ey.etf=etf;vector.Ey.ep=ep;
vector.Ey.ex0=ex0;vector.Ey.eu0=eu0;
vector.Ey.exf=exf;vector.Ey.euf=euf;


% dfdz where z=[tf p x0 ...u0... xf uf];
% -------------------------------------
Fxu=[kron(speye((M-1)/N),dfdx),repmat(dfdu,(M-1)/N,1)];
df=[repmat(ones(n,nt),M-1,1) repmat(dfdp,M-1,1) kron(speye(N),Fxu)];
df=[df zeros((M-1)*n,n); ones(n,nt) dfdp zeros(n,size(df,2)-m-nt-np) dfdu dfdx];

colFxu=[repmat(overlapping(Fxu),N,1);overlapping(dfdx)];
if isempty(dfdp); dp=0;else dp=max(overlapping(dfdp));end


colgroup=[ones(nt);ones(~(np==0))*(overlapping(dfdp)+nt);...      
         nt+dp+colFxu];ngcol=max(colgroup);                      


eg=spalloc(nz,ngcol,nz);
vector.f.etf=zeros(1,ngcol);vector.f.etf(1)=nt;
vector.f.ep=eg;vector.f.ex=eg;vector.f.eu=eg;

ixf=spalloc(n*M,ngcol,n*M*ngcol);


for i=1:ngcol

eg(colgroup==i,i)=1;
 
    
[a,b]=find(df*spdiags(eg(:,i),0,nz,nz));
ixf(a,i)=b;
ixf(ixf(:,i)==0,i)=1;
end


ex=data.map.Vx*eg;eu=data.map.Vu*eg;


vector.f.ex=cell(ngcol,1);vector.f.eu=cell(ngcol,1);vector.f.ep=cell(ngcol,1);


for i=1:ngcol
   
    vector.f.ex{i}=reshape(ex(:,i)',n,M)';
    eus=kron(reshape(eu(:,i)',m,N)',ones((M-1)/N,1)); 
    vector.f.eu{i}=[eus;eus(end,:)];
    vector.f.ep{i}=repmat(eg(nt+1:nt+np,i)',M,1);
end

index.f=ixf;



% dgdz where z=[tf p x0 ...u0...  uf xf]; 

if ng
Gxu=[kron(speye((M-1)/N),sparse(sparsity.dgdx)),repmat(sparse(sparsity.dgdu),(M-1)/N,1)];

dg=[repmat(dgdt,M,1) repmat(dgdp,M,1) [kron(speye(N),Gxu),  zeros(ng*(M-1),n);sparsity.dgdt sparsity.dgdp zeros(ng,(M-1)*n+(N-1)*m-nt-np),  sparsity.dgdu sparsity.dgdx]];

   
cGxu=overlapping(Gxu);   
colGxu=[repmat(cGxu,N,1); cGxu(1:n)];
if isempty(dgdp); dp=0;else dp=max(overlapping(dgdp));end


colgroup=[ones(nt);ones(~(np==0))*(overlapping(dgdp)+nt);...
         nt+dp+colGxu];ngcol=max(colgroup);



eg=spalloc(nz,ngcol,ngcol*nz);
vector.g.etf=zeros(1,ngcol);vector.g.etf(1)=nt;
vector.g.ep=eg;vector.g.ex=eg;vector.g.eu=eg;

ig=spalloc(ng*M,ngcol,ng*M*ngcol);


for i=1:ngcol
    eg(colgroup==i,i)=1;
    [a,b]=find(dg*spdiags(eg(:,i),0,nz,nz));
    
    ig(a,i)=b;
       
    ig(ig(:,i)==0,i)=1;
  
end

ex=data.map.Vx*eg;eu=data.map.Vu*eg;

vector.g.ex=cell(ngcol,1);vector.g.eu=cell(ngcol,1);vector.g.ep=cell(ngcol,1);

for i=1:ngcol
    vector.g.ex{i}=reshape(ex(:,i)',n,M)';
    ugs=kron(reshape(eu(:,i)',m,N)',ones((M-1)/N,1));
    vector.g.eu{i}=[ugs;ugs(end,:)]; % can be optimized
    vector.g.ep{i}=repmat(eg(nt+1:nt+np,i)',M,1);
end

index.g=ig;
end




% dbdzb where zb=[tf p x0 u0 uf xf];
% -----------------------------------



db=[dbdtf dbdp dbdx0 dbdu0 dbduf dbdxf];
colgroup=overlapping(db);nbcol=max(colgroup);



eb=zeros(2*(n+m)+np+nt,nbcol);
vector.b.etf=eb(1,:);vector.b.ep=eb(1:np,:);
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

if nt; vector.b.etf=eb(nt,:); end
if np; vector.b.ep=eb(nt+1:nt+np,:); end

vector.b.ex0=eb(nt+np+1:nt+np+n,:);
vector.b.eu0=eb(nt+np+n+1:nt+np+n+m,:);
vector.b.euf=eb(nt+np+n+m+1:nt+np+n+2*m,:);
vector.b.exf=eb(nt+np+n+2*m+1:nt+np+2*n+2*m,:);

end


% shift index
ibm=ib;
ibm((ib>nt+np+n)&(ib<nt+np+n+m+1))=ib((ib>nt+np+n)&(ib<nt+np+n+m+1))+((M-1)/N-1)*n;   
ibm(ib>nt+np+n+m)=nz-nt-np-2*n-2*m+ib((ib>nt+np+n+m));
index.idb=ib;
index.b=ibm;







%------------- END OF CODE --------------
  