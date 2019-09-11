function [HL_new,HE_new,Hf_new,Hg_new,Hb_new] = convertHessianANStruct(HL,HE,Hf,Hg,Hb,data)
%convertHessianANStruct - convert the Hessian structure to a format compatible with h method formulation of ICLOCS2 
%
% Syntax:  [HL_new,HE_new,Hf_new,Hg_new,Hb_new] = convertHessianANStruct(HL,HE,Hf,Hg,Hb,data)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
% 
%------------- BEGIN CODE --------------

[nt,np,n,m,ng,nb]=deal(data.sizes{1:6});

h_method=~(strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR'));

if ~isempty(HL) && h_method
    HL_new=cell(size(HL));
    HL_new(1:nt,1:nt)=HL(end-nt+1:end,end-nt+1:end);
    HL_new(1+nt:nt+np,1+nt:nt+np)=HL(end-nt-np+1:end-nt,end-nt-np+1:end-nt);
    HL_new(1+nt+np:end,1+nt+np:end)=HL(1:end-nt-np,1:end-nt-np);
    HL_new(1:nt,1+nt:nt+np)=transpose(HL(end-nt-np+1:end-nt,end-nt+1:end));
    HL_new(1:nt,1+nt+np:end)=transpose(HL(1:end-nt-np,end-nt+1:end));
    HL_new(1+nt:nt+np,1+nt+np:end)=transpose(HL(1:end-nt-np,end-nt-np+1:end-nt));
else
    HL_new=HL;
end

if ~isempty(HE) && h_method
    HE_new=cell(size(HE));
    HE_new(1:nt,1:nt)=HE(end-nt+1:end,end-nt+1:end);%t t
    HE_new(1+nt:nt+np,1+nt:nt+np)=HE(end-nt-np+1:end-nt,end-nt-np+1:end-nt);%p p
    HE_new(1+nt+np:nt+np+n,1+nt+np:nt+np+n)=HE(1:n,1:n); %x0 x0
    HE_new(1+nt+np+n:nt+np+n+m,1+nt+np+n:nt+np+n+m)=HE(1+n*2:n*2+m,1+n*2:n*2+m); %u0 u0
    HE_new(1+nt+np+n+m:nt+np+n*2+m,1+nt+np+n+m:nt+np+n*2+m)=HE(1+n:n*2,1+n:n*2); %xf xf
    HE_new(1+nt+np+n*2+m:nt+np+n*2+m*2,1+nt+np+n*2+m:nt+np+n*2+m*2)=HE(1+n*2+m:n*2+m*2,1+n*2+m:n*2+m*2); %uf uf
    
    HE_new(1:nt,1+nt:nt+np)=transpose(HE(end-nt-np+1:end-nt,end-nt+1:end)); %t p
    HE_new(1:nt,1+nt+np:nt+np+n)=transpose(HE(1:n,end-nt+1:end)); %t x0
    HE_new(1:nt,1+nt+np+n:nt+np+n+m)=transpose(HE(1+n*2:n*2+m,end-nt+1:end)); %t u0
    HE_new(1:nt,1+nt+np+n+m:nt+np+n*2+m)=transpose(HE(1+n:n*2,end-nt+1:end)); %t xf
    HE_new(1:nt,1+nt+np+n*2+m:nt+np+n*2+m*2)=transpose(HE(1+n*2+m:n*2+m*2,end-nt+1:end)); %t uf
    
    HE_new(1+nt:nt+np,1+nt+np:nt+np+n)=transpose(HE(1:n,end-nt-np+1:end-nt)); %p x0
    HE_new(1+nt:nt+np,1+nt+np+n:nt+np+n+m)=transpose(HE(1+n*2:n*2+m,end-nt-np+1:end-nt)); %p u0
    HE_new(1+nt:nt+np,1+nt+np+n+m:nt+np+n*2+m)=transpose(HE(1+n:n*2,end-nt-np+1:end-nt)); %p xf
    HE_new(1+nt:nt+np,1+nt+np+n*2+m:nt+np+n*2+m*2)=transpose(HE(1+n*2+m:n*2+m*2,end-nt-np+1:end-nt)); %p uf
    
    HE_new(1+nt+np:nt+np+n,1+nt+np+n:nt+np+n+m)=HE(1:n,1+n*2:n*2+m); %x0 u0
    HE_new(1+nt+np:nt+np+n,1+nt+np+n+m:nt+np+n*2+m)=HE(1:n,1+n:n*2); %x0 xf
    HE_new(1+nt+np:nt+np+n,1+nt+np+n*2+m:nt+np+n*2+m*2)=HE(1:n,1+n*2+m:n*2+m*2); %x0 uf
    
    HE_new(1+nt+np+n:nt+np+n+m,1+nt+np+n+m:nt+np+n*2+m)=transpose(HE(1+n:n*2,1+n*2:n*2+m)); %u0 xf
    HE_new(1+nt+np+n:nt+np+n+m,1+nt+np+n*2+m:nt+np+n*2+m*2)=HE(1+n*2:n*2+m,1+n*2+m:n*2+m*2); %u0 uf
    
    HE_new(1+nt+np+n+m:nt+np+n*2+m,1+nt+np+n*2+m:nt+np+n*2+m*2)=transpose(HE(1+n:n*2,1+n*2+m:n*2+m*2)); %xf uf
else
    HE_new=HE;
end

if ~isempty(Hf) && h_method
    Hf_new=cell(size(Hf));
    Hf_new(1:nt,1:nt)=Hf(end-nt+1:end,end-nt+1:end);
    Hf_new(1+nt:nt+np,1+nt:nt+np)=Hf(end-nt-np+1:end-nt,end-nt-np+1:end-nt);
    Hf_new(1+nt+np:end,1+nt+np:end)=Hf(1:end-nt-np,1:end-nt-np);
    Hf_new(1:nt,1+nt:nt+np)=transpose(Hf(end-nt-np+1:end-nt,end-nt+1:end));
    Hf_new(1:nt,1+nt+np:end)=transpose(Hf(1:end-nt-np,end-nt+1:end));
    Hf_new(1+nt:nt+np,1+nt+np:end)=transpose(Hf(1:end-nt-np,end-nt-np+1:end-nt));
else
    Hf_new=Hf;
end

if ~isempty(Hg) && h_method
    Hg_new=cell(size(Hg));
    Hg_new(1:nt,1:nt)=Hg(end-nt+1:end,end-nt+1:end);
    Hg_new(1+nt:nt+np,1+nt:nt+np)=Hg(end-nt-np+1:end-nt,end-nt-np+1:end-nt);
    Hg_new(1+nt+np:end,1+nt+np:end)=Hg(1:end-nt-np,1:end-nt-np);
    Hg_new(1:nt,1+nt:nt+np)=transpose(Hg(end-nt-np+1:end-nt,end-nt+1:end));
    Hg_new(1:nt,1+nt+np:end)=transpose(Hg(1:end-nt-np,end-nt+1:end));
    Hg_new(1+nt:nt+np,1+nt+np:end)=transpose(Hg(1:end-nt-np,end-nt-np+1:end-nt));
else
    Hg_new=Hg;
end

if ~isempty(Hb) && h_method
    Hb_new=cell(size(Hb));
    Hb_new(1:nt,1:nt)=Hb(end-nt+1:end,end-nt+1:end);%t t
    Hb_new(1+nt:nt+np,1+nt:nt+np)=Hb(end-nt-np+1:end-nt,end-nt-np+1:end-nt);%p p
    Hb_new(1+nt+np:nt+np+n,1+nt+np:nt+np+n)=Hb(1:n,1:n); %x0 x0
    Hb_new(1+nt+np+n:nt+np+n+m,1+nt+np+n:nt+np+n+m)=Hb(1+n*2:n*2+m,1+n*2:n*2+m); %u0 u0
    Hb_new(1+nt+np+n+m:nt+np+n*2+m,1+nt+np+n+m:nt+np+n*2+m)=Hb(1+n:n*2,1+n:n*2); %xf xf
    Hb_new(1+nt+np+n*2+m:nt+np+n*2+m*2,1+nt+np+n*2+m:nt+np+n*2+m*2)=Hb(1+n*2+m:n*2+m*2,1+n*2+m:n*2+m*2); %uf uf
    
    Hb_new(1:nt,1+nt:nt+np)=transpose(Hb(end-nt-np+1:end-nt,end-nt+1:end)); %t p
    Hb_new(1:nt,1+nt+np:nt+np+n)=transpose(Hb(1:n,end-nt+1:end)); %t x0
    Hb_new(1:nt,1+nt+np+n:nt+np+n+m)=transpose(Hb(1+n*2:n*2+m,end-nt+1:end)); %t u0
    Hb_new(1:nt,1+nt+np+n+m:nt+np+n*2+m)=transpose(Hb(1+n:n*2,end-nt+1:end)); %t xf
    Hb_new(1:nt,1+nt+np+n*2+m:nt+np+n*2+m*2)=transpose(Hb(1+n*2+m:n*2+m*2,end-nt+1:end)); %t uf
    
    Hb_new(1+nt:nt+np,1+nt+np:nt+np+n)=transpose(Hb(1:n,end-nt-np+1:end-nt)); %p x0
    Hb_new(1+nt:nt+np,1+nt+np+n:nt+np+n+m)=transpose(Hb(1+n*2:n*2+m,end-nt-np+1:end-nt)); %p u0
    Hb_new(1+nt:nt+np,1+nt+np+n+m:nt+np+n*2+m)=transpose(Hb(1+n:n*2,end-nt-np+1:end-nt)); %p xf
    Hb_new(1+nt:nt+np,1+nt+np+n*2+m:nt+np+n*2+m*2)=transpose(Hb(1+n*2+m:n*2+m*2,end-nt-np+1:end-nt)); %p uf
    
    Hb_new(1+nt+np:nt+np+n,1+nt+np+n:nt+np+n+m)=Hb(1:n,1+n*2:n*2+m); %x0 u0
    Hb_new(1+nt+np:nt+np+n,1+nt+np+n+m:nt+np+n*2+m)=Hb(1:n,1+n:n*2); %x0 xf
    Hb_new(1+nt+np:nt+np+n,1+nt+np+n*2+m:nt+np+n*2+m*2)=Hb(1:n,1+n*2+m:n*2+m*2); %x0 uf
    
    Hb_new(1+nt+np+n:nt+np+n+m,1+nt+np+n+m:nt+np+n*2+m)=transpose(Hb(1+n:n*2,1+n*2:n*2+m)); %u0 xf
    Hb_new(1+nt+np+n:nt+np+n+m,1+nt+np+n*2+m:nt+np+n*2+m*2)=Hb(1+n*2:n*2+m,1+n*2+m:n*2+m*2); %u0 uf
    
    Hb_new(1+nt+np+n+m:nt+np+n*2+m,1+nt+np+n*2+m:nt+np+n*2+m*2)=transpose(Hb(1+n:n*2,1+n*2+m:n*2+m*2)); %xf uf
else
    Hb_new=Hb;
end

end

