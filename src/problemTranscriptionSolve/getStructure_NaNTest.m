function [ sparsity ] = getStructure_NaNTest( problem, options, data, sparsity_num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

sparsity=sparsity_num;

if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
    [nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:17});
else
    [nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:13});
end
nrc=nrcl+nrcu+nrce;

[L,E,f,g,avrc,b]=deal(data.functions{:});

X_test=rand(1,n)*0.5;
U_test=rand(1,m)*0.5;
P_test=rand(1,np)*0.5;
X_test(X_test<problem.states.xl)=problem.states.xl(X_test<problem.states.xl);
X_test(X_test>problem.states.xu)=problem.states.xu(X_test>problem.states.xu);
U_test(U_test<problem.inputs.ul)=problem.inputs.ul(U_test<problem.inputs.ul);
U_test(U_test>problem.inputs.uu)=problem.inputs.uu(U_test>problem.inputs.uu);
P_test(P_test<problem.parameters.pl)=problem.parameters.pl(P_test<problem.parameters.pl);
P_test(P_test>problem.parameters.pu)=problem.parameters.pu(P_test>problem.parameters.pu);

dfdx=zeros(n,n);
dgdx=zeros(ng,n);
dLdx=zeros(1,n);
dEdx0=zeros(1,n);
dEdxf=zeros(1,n);
for i=1:n
    X_temp=X_test;
    X_temp(1,i)=NaN;
    dfdx(:,i)=f(X_temp,U_test,P_test,1,problem.data)';
    dLdx(1,i)=L(X_temp,[],U_test,[],P_test,1,problem.data);
    dEdx0(1,i)=E(X_temp,X_test,U_test,U_test,P_test,1,1,problem.data);
    dEdxf(1,i)=E(X_test,X_temp,U_test,U_test,P_test,1,1,problem.data);
    if ng
        dgdx(:,i)=g(X_temp,U_test,P_test,1,problem.data)';
    end
end
dfdx(~isnan(dfdx))=0;
dgdx(~isnan(dgdx))=0;
dLdx(~isnan(dLdx))=0;
dEdx0(~isnan(dEdx0))=0;
dEdxf(~isnan(dEdxf))=0;
dfdx(isnan(dfdx))=1;
dgdx(isnan(dgdx))=1;
dLdx(isnan(dLdx))=1;
dEdx0(isnan(dEdx0))=1;
dEdxf(isnan(dEdxf))=1;

dfdu=zeros(n,m);
dgdu=zeros(ng,m);
dLdu=zeros(1,m);
dEdu0=zeros(1,m);
dEduf=zeros(1,m);
for i=1:m
    U_temp=U_test;
    U_temp(1,i)=NaN;
    dfdu(:,i)=f(X_test,U_temp,P_test,1,problem.data)';
    dLdu(1,i)=L(X_test,[],U_temp,[],P_test,1,problem.data);
    dEdu0(1,i)=E(X_test,X_test,U_temp,U_test,P_test,1,1,problem.data);
    dEduf(1,i)=E(X_test,X_test,U_test,U_temp,P_test,1,1,problem.data);
    if ng
        dgdu(:,i)=g(X_test,U_temp,P_test,1,problem.data)';
    end
end
dfdu(~isnan(dfdu))=0;
dgdu(~isnan(dgdu))=0;
dLdu(~isnan(dLdu))=0;
dEdu0(~isnan(dEdu0))=0;
dEduf(~isnan(dEduf))=0;
dfdu(isnan(dfdu))=1;
dgdu(isnan(dgdu))=1;
dLdu(isnan(dLdu))=1;
dEdu0(isnan(dEdu0))=1;
dEduf(isnan(dEduf))=1;

sparsity.dfdx=sparse(dfdx);
sparsity.dgdx=sparse(dgdx);
sparsity.dLdx=sparse(dLdx);
sparsity.dEdx0=sparse(dEdx0);
sparsity.dEdxf=sparse(dEdxf);

sparsity.dfdu=sparse(dfdu);
sparsity.dgdu=sparse(dgdu);
sparsity.dLdu=sparse(dLdu);
sparsity.dEdu0=sparse(dEdu0);
sparsity.dEduf=sparse(dEduf);


end

