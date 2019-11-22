function [jac,Jf]=jacConstz(df,dg,g,f,avrc,X,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data)
%JACCONSTZ - Return the Jacobian of the Constraints with respect to the 
%            variable z  when the analytic option has been selected
%
% Syntax:  [jac,Jf]=jacConstz(df,dg,g,f,avrc,X,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data)
%
% Inputs:
%    see directCollocation.m
%
% Outputs:
%    jac, Jf - jacobian information
%
% Other m-files required: none
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

%------------- BEGIN CODE ---------------


e=data.options.perturbation.J;                                % pertubation size
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,~]=deal(data.sizes{1:13});
nrc=nrcl+nrcu+nrce;
nz=nt+np+M*n+N*m;
vdat=data.data;
% t0=data.t0;
% k0=data.k0;


% Compute fz and gz
%------------

fz=spalloc(n*M,nz,n*M*(n+m)+np+nt);
gz=spalloc(ng*M,nz,ng*M*(n+m+np+nt));
Jaf=[];

if df.flag && ng && dg.flag
    [ fz,Jf ] = jacConst_AN_F( df, fz, M, n, m, np, nt, nz, f, X, U, P, t0, tf, T, e, vdat, data );
    [ gz ] = jacConst_AN_G( gz, M, nt, np, ng, nz, T, data );
elseif ~df.flag && ng && ~dg.flag
    [ fz,gz,Jf ] = jacConst_FD_FG( fz, gz, M, n, ng, nz, fg, X, U, P, t0, tf, T, e, vdat, data );
elseif df.flag 
    [ fz,Jf ] = jacConst_AN_F( df, fz, M, n, m, np, nt, nz, f, X, U, P, t0, tf, T, e, vdat, data );
    if ng
        if dg.flag
            [ gz ] = jacConst_AN_G( gz, M, nt, np, ng, nz, T, data );
        else
            [ gz ] = jacConst_FD_G( gz, M, ng, nz, g, X, U, P, t0, tf, T, e, vdat, data );
        end
    end
elseif ~df.flag
    [ fz,Jf ] = jacConst_FD_F( fz, M, n, nz, f, X, U, P, t0, tf, T, e, vdat, data );
    if ng
        if dg.flag
            [ gz ] = jacConst_AN_G( gz, M, nt, np, ng, nz, T, data );
        else
            [ gz ] = jacConst_FD_G( gz, M, ng, nz, g, X, U, P, t0, tf, T, e, vdat, data );
        end
    end
end

% Compute rcz
%------------
rcz=spalloc(nrc,nz,nrc*(n+m+np+nt));

if nrc
    [ rcz ] = jacConst_RC( rcz, nrc, nz, avrc, X, U, P, t0, tf, T, e, data );
end

% Compute bz
%------------

bz=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
if nb

    if db.flag 
        [ bz ] = jacConst_AN_B( bz, db, nt, data );
    else     
        [ bz ] = jacConst_FD_B( bz, nb, nz, b, x0, xf, u0, uf, p, t0, tf, e, vdat, data );
    end
end

% Map derivatives to the jacobian
%---------------------------------
jac=[[sparse(n,nt) sparse(n,np) speye(n), sparse(n,(M-1)*n+N*m)]*data.cx0;data.map.A*data.map.Vx+data.map.B*fz;gz;rcz;bz];






%------------- END CODE --------------