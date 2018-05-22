function [ fx, fxx ] = getAdigator4ICLOCS( X, U, t, p, probinfo )
%genAdigator4ICLOCS - Initialize Adigator for use with ICLOCS
%
% Syntax:  [ ] = genAdigator4ICLOCS( options, probinfo, n, m, nt, np )
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    probinfo - structured variable containing the values of additional data used inside
%          the function 
%
% Outputs:
%    fx, fxx - First and second derivative infomation
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

n=probinfo.sizes{3};
m=probinfo.sizes{4};

if probinfo.options.scaling
    X=scale_variables_back( X, probinfo.data.Xscale, probinfo.data.Xshift );
    U=scale_variables_back( U, probinfo.data.Uscale, probinfo.data.Ushift );
end

Xvec.f = X;
Xvec.dY = ones(size(Xvec.f));
Uvec.f  = U;
Uvec.dY = ones(size(Uvec.f));
Tvec.f = t;
Tvec.dT = ones(size(Tvec.f));
Pvec.f = p;

% Get Vectorized Derivs
if strcmp(probinfo.options.ipopt.hessian_approximation,'limited-memory')
    Fvec = dynamics_Y(Xvec,Uvec,Pvec,Tvec,probinfo.data);
else
    Fvec = dynamics_YY(Xvec,Uvec,Pvec,Tvec,probinfo.data);
end

if length(Fvec.dY_size)==1
    Fvec_dY_size1=size(Fvec.dY,1);
    Fvec_dY_size2=n+m;
    Fvec_dY_location1=1;
    Fvec_dY_location2=Fvec.dY_location(:,1);
else
    Fvec_dY_size1=size(Fvec.dY,1);
    Fvec_dY_size2=Fvec.dY_size(2);
    Fvec_dY_location1=Fvec.dY_location(:,1);
    Fvec_dY_location2=Fvec.dY_location(:,2);
end
if length(Fvec.dYdY_size)==2
    Fvec_dYdY_size3=n+m;
    Fvec_dYdY_size=[1,Fvec.dYdY_size(1),Fvec_dYdY_size3];
    Fvec_dYdY_location3=Fvec.dYdY_location(:,2);
else
    Fvec_dYdY_size3=Fvec.dYdY_size(3);
    Fvec_dYdY_size=Fvec.dYdY_size;
    Fvec_dYdY_location3=Fvec.dYdY_location(:,3);
end

Fx = cell(1,Fvec_dY_size1);fx = cell(1,Fvec_dY_size2);
if strcmp(probinfo.options.ipopt.hessian_approximation,'limited-memory') || ~isfield(Fvec,'dYdY')
    for i=1:Fvec_dY_size1
        Fx{i}=full(sparse(Fvec.dY_location(:,1),Fvec_dY_location2,Fvec.dY(i,:)));
    end
    for i=1:Fvec_dY_size1
        for j=1:Fvec_dY_size2
            fx{j}(:,i)=Fx{i}(:,j);
        end
    end
    fxx=[];
else
    Fxx = cell(1,size(Fvec.dYdY,1));Fxx_2d = cell(size(Fvec.dYdY,1),Fvec_dYdY_size3); fxx = cell(Fvec.dYdY_size(2),Fvec_dYdY_size3);
    for i=1:Fvec_dY_size1
        Fx{i}=full(sparse(Fvec_dY_location1,Fvec_dY_location2,Fvec.dY(i,:)));
        Fxx3dmat=zeros(Fvec_dYdY_size);
        idx = sub2ind(size(Fxx3dmat), Fvec.dYdY_location(:,1),Fvec.dYdY_location(:,2),Fvec_dYdY_location3);
        Fxx3dmat(idx)=Fvec.dYdY(i,:);
        Fxx{i}=Fxx3dmat;
    end
    for i=1:Fvec_dY_size1
        for j=1:Fvec_dY_size2
            fx{j}(:,i)=Fx{i}(:,j);
        end
    end
    for i=1:size(Fvec.dYdY,1)
        for j=1:Fvec_dYdY_size3
            Fxx_2d{i,j}=Fxx{1,i}(:,:,j);
        end
    end
    for k=1:Fvec_dYdY_size3
        for i=1:size(Fvec.dYdY,1)
            for j=1:Fvec.dYdY_size(2)
                fxx{j,k}(:,i)=Fxx_2d{i,k}(:,j);
            end
        end
    end
end

if probinfo.options.scaling
    for i=1:Fvec_dY_size2
        for j=1:Fvec.dY_size(1)
            if i<=n
                fx{i}(j,:)=fx{i}(j,:)/(probinfo.data.Xscale(i)/probinfo.data.Xscale(j));
            elseif i<=(m+n)
                fx{i}(j,:)=fx{i}(j,:)/(probinfo.data.Uscale(i-n)/probinfo.data.Xscale(j));
            end
        end
    end
    if ~strcmp(probinfo.options.ipopt.hessian_approximation,'limited-memory') && isfield(Fvec,'dYdY')
        for k=1:Fvec_dYdY_size3
            for i=1:Fvec.dYdY_size(2)
                for j=1:Fvec.dY_size(1)
                    if i<=n && k<=n
                        fxx{i,k}(j,:)=fxx{i,k}(j,:)/(probinfo.data.Xscale(i)*probinfo.data.Xscale(k)/probinfo.data.Xscale(j));
                    elseif i<=(m+n) && k<=n
                        fxx{i,k}(j,:)=fxx{i,k}(j,:)/(probinfo.data.Uscale(i-n)*probinfo.data.Xscale(k)/probinfo.data.Xscale(j));
                    elseif i<=n && k<=(m+n)
                        fxx{i,k}(j,:)=fxx{i,k}(j,:)/(probinfo.data.Xscale(i)*probinfo.data.Uscale(k-n)/probinfo.data.Xscale(j));
                    elseif i<=(m+n) && k<=(m+n)
                        fxx{i,k}(j,:)=fxx{i,k}(j,:)/(probinfo.data.Uscale(i-n)*probinfo.data.Uscale(k-n)/probinfo.data.Xscale(j));
                    end
                end
            end
        end
    end
end


end

