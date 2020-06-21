function [jac,Jf]=jacConstzLGR(df,dg,g,f,X_Np1,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data)

%jacConstzLGR - Return the Jacobian of the Constraints with respect to the 
%            variable z  when the analytic option has been selected
%
% Syntax:  [jac,Jf]=jacConstzLGR(df,dg,g,f,X,U,P,T,db,b,x0,xf,u0,uf,p,t0,tf,data)
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
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,~]=deal(data.sizes{1:17});
nrc=nrcl+nrcu+nrce;
nz=nt+np+(M+1)*n+M*m;
vdat=data.data;
fg=vdat.functionfg;
X=X_Np1(1:M,:);
DTLP=repmat(data.t_segment_end,1,n);
DT_ratio_diff=repmat(data.tau_segment_ratio_diff,1,n);

% Compute fz and gz
%------------

fz=spalloc(n*M,nz,data.map.spmatsize.jSf);
gz=spalloc(ng*M,nz,data.map.spmatsize.jSg);
Jaf=[];

if df.flag && ng && dg.flag
    [ fz,Jf ] = jacConst_LGR_AN_F( df, fz, M, n, m, np, nt, nz, f, X, U, P, t0, tf, T, e, DT_ratio_diff, DTLP, vdat, data );
    [ gz ] = jacConst_LGR_AN_G( dg, gz, M, nt, np, ng, nz, T, data );
elseif ~df.flag && ng && ~dg.flag && size(data.FD.index.f,2)==size(data.FD.index.g,2)
    [ fz,gz,Jf ] = jacConst_LGR_FD_FG( fz, gz, M, n, ng, nz, fg, X, U, P, t0, tf, T, e, DT_ratio_diff, DTLP, vdat, data );
elseif df.flag 
    [ fz,Jf ] = jacConst_LGR_AN_F( df, fz, M, n, m, np, nt, nz, f, X, U, P, t0, tf, T, e, DT_ratio_diff, DTLP, vdat, data );
    if ng
        if dg.flag
            [ gz ] = jacConst_LGR_AN_G( dg, gz, M, nt, np, ng, nz, T, data );
        else
            [ gz ] = jacConst_LGR_FD_G( gz, M, ng, nz, g, X, U, P, t0, tf, T, e, vdat, data );
        end
    end
elseif ~df.flag
    [ fz,Jf ] = jacConst_LGR_FD_F( fz, M, n, nz, f, X, U, P, t0, tf, T, e, DT_ratio_diff, DTLP, vdat, data );
    if ng
        if dg.flag
            [ gz ] = jacConst_LGR_AN_G( dg, gz, M, nt, np, ng, nz, T, data );
        else
            [ gz ] = jacConst_LGR_FD_G( gz, M, ng, nz, g, X, U, P, t0, tf, T, e, vdat, data );
        end
    end
end


% Compute rcz
%------------

rcz=spalloc(nrc,nz,data.map.spmatsize.jSrc);
if nrc
    [ rcz ] = jacConst_LGR_RC( rcz, nrc, nz, avrc, X_Np1, U, P, t0, tf, T, e, data );
end


% Compute bz
%------------

bz=spalloc(nb,nz,data.map.spmatsize.jSb);
if nb

    if db.flag 
        [ bz ] = jacConst_LGR_AN_B( bz, db, nt, data );
    else     
        [ bz ] = jacConst_LGR_FD_B( bz, nb, nz, b, x0, xf, u0, uf, p, t0, tf, e, vdat, data );
    end
end


% Map derivatives to the jacobian
%---------------------------------
jac=[[kron(speye(n),data.map.D_structure) zeros(M*n,M*m+np+nt) ]+fz;gz(data.gAllidx,:);rcz;bz];






%------------- END CODE --------------