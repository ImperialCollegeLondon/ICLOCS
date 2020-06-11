function jacConst_mp=constraintJacobian_mp(z_mp,auxdata)

%CONSTRAINTJACOBIAN - Evaluate the Jacobian of the constraints, for a
%multiphase problem
%
% Syntax: jacConst_mp=constraintJacobian_mp(z_mp,auxdata)
%
% Inputs:
%    z - NLP variable 
%    data - Structure with data needed to evaluate CostFn
%
% Outputs:
%    jacConst - Sparse matrix with jacobian of the constraints at z
%
% Other m-files required: multipleShooting.m or directTranscription.m
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


%------------- BEGIN CODE --------------
nbl=auxdata.mpdata.mpsizes.nbl;

jacConst_mp=sparse(auxdata.mpdata.mpsizes.nConst-nbl,auxdata.mpdata.mpsizes.nz);

x0=cell(auxdata.mpdata.mpsizes.nphase,1);xf=cell(auxdata.mpdata.mpsizes.nphase,1);
u0=cell(auxdata.mpdata.mpsizes.nphase,1);uf=cell(auxdata.mpdata.mpsizes.nphase,1);

for i=1:length(auxdata.phasedata)
    data=auxdata.phasedata{i};
    z=zeros(data.nz,1);
    z(data.zidx.org.z,1)=z_mp(data.zidx.mp.z,1);

    switch data.options.transcription

        case {'multiple_shooting'}

             jacConst=multipleShooting('jacConst',z,data);

        case {'direct_collocation'}

            switch data.options.discretization

               case{'globalLGR','hpLGR'}

                     [jacConst,XU0f]=directCollocationLGR('jacConst',z,data,i);
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;
                     
                case {'resMinConstructionForSolution'}

                    jacConst=resMinConstructionForSolution('jacConst',z,data);

                case {'resMinInterpolationForSolution'}

    %                 jacConst=resMinInterpolationForSolution('jacConst',z,data);
                    [jacConst,XU0f]=resMinOptimization('jacConst',z,data);
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;

                otherwise % Direct Transcription Method

                     [jacConst,XU0f]=directCollocation('jacConst',z,data,i);
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;
            end

        case {'integral_res_min'}

            [jacConst,XU0f]=resMinOptimization('jacConst',z,data);
             x0{i}=XU0f.x0;
             xf{i}=XU0f.xf;
             u0{i}=XU0f.u0;
             uf{i}=XU0f.uf;
    end
    
%     jacConst_mp(data.jSidx.mp.all.row,data.jSidx.mp.all.col)=jacConst_mp(data.jSidx.mp.all.row,data.jSidx.mp.all.col)+jacConst(data.jSidx.org.all.row,data.jSidx.org.all.col);
    jacConst_mp(data.jSidx.mp.XUg.row,data.jSidx.mp.XUg.col)=jacConst_mp(data.jSidx.mp.XUg.row,data.jSidx.mp.XUg.col)+jacConst(data.jSidx.org.XUg.row,data.jSidx.org.XUg.col);
    jacConst_mp(data.jSidx.mp.P.row,data.jSidx.mp.P.col)=jacConst_mp(data.jSidx.mp.P.row,data.jSidx.mp.P.col)+jacConst(data.jSidx.org.P.row,data.jSidx.org.P.col);
    jacConst_mp(data.jSidx.mp.T.row,data.jSidx.mp.T.col)=jacConst_mp(data.jSidx.mp.T.row,data.jSidx.mp.T.col)+jacConst(data.jSidx.org.T.row,data.jSidx.org.T.col);
    jacConst_mp(data.jSidx.mp.B_XUg.row,data.jSidx.mp.B_XUg.col)=jacConst_mp(data.jSidx.mp.B_XUg.row,data.jSidx.mp.B_XUg.col)+jacConst(data.jSidx.org.B_XUg.row,data.jSidx.org.B_XUg.col);
    jacConst_mp(data.jSidx.mp.B_P.row,data.jSidx.mp.B_P.col)=jacConst_mp(data.jSidx.mp.B_P.row,data.jSidx.mp.B_P.col)+jacConst(data.jSidx.org.B_P.row,data.jSidx.org.B_P.col);
    jacConst_mp(data.jSidx.mp.B_T.row,data.jSidx.mp.B_T.col)=jacConst_mp(data.jSidx.mp.B_T.row,data.jSidx.mp.B_T.col)+jacConst(data.jSidx.org.B_T.row,data.jSidx.org.B_T.col);
end

linkfunctions=auxdata.mpdata.linkfunctions;
mpdata=auxdata.mpdata;
t0=z_mp(end-mpdata.mpsizes.nt+1:end-1);
tf=z_mp(end-mpdata.mpsizes.nt+2:end);
p=z_mp(end-mpdata.mpsizes.nt-mpdata.mpsizes.np+1:end-mpdata.mpsizes.nt);
if isfield(data.data,'Pscale')
   p=scale_variables_back( p', data.data.Pscale, data.data.Pshift )';
end

blz=spalloc(nbl,auxdata.mpdata.mpsizes.nz,(auxdata.mpdata.mpsizes.nxu0f_inc(end)+auxdata.mpdata.mpsizes.nt+auxdata.mpdata.mpsizes.np)*nbl);
if nbl
    idx=auxdata.mpdata.linkConst.all.idx;
    nfd=size(idx,2);
    % ez=auxdata.mpdata.options.mp.perturbation.J*auxdata.mpdata.map.all;
    for i=1:nfd
        x0p=x0;x0m=x0;
        xfp=xf;xfm=xf;
        u0p=u0;u0m=u0;
        ufp=uf;ufm=uf;

        if mpdata.map.allidx(i)~=0
            x0p{mpdata.map.allidx(i)}=x0p{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.ex0(:,i);
            x0m{mpdata.map.allidx(i)}=x0m{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.ex0(:,i);
            xfp{mpdata.map.allidx(i)}=xfp{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.exf(:,i);
            xfm{mpdata.map.allidx(i)}=xfm{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.exf(:,i);
            u0p{mpdata.map.allidx(i)}=u0p{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.eu0(:,i);
            u0m{mpdata.map.allidx(i)}=u0m{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.eu0(:,i);
            ufp{mpdata.map.allidx(i)}=ufp{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.euf(:,i);
            ufm{mpdata.map.allidx(i)}=ufm{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.euf(:,i);
        end

        [blp_l, blp_nl]=linkfunctions(x0p,xfp,u0p,ufp,p+mpdata.map.ep(:,i),t0+mpdata.map.et0(:,i),tf+mpdata.map.etf(:,i),mpdata.data);
        [blm_l, blm_nl]=linkfunctions(x0m,xfm,u0m,ufm,p-mpdata.map.ep(:,i),t0-mpdata.map.et0(:,i),tf-mpdata.map.etf(:,i),mpdata.data);
        blz=blz+sparse(1:nbl,idx(:,i),([blp_l;blp_nl]-[blm_l;blm_nl])/(2*auxdata.mpdata.options.mp.perturbation.J),nbl,auxdata.mpdata.mpsizes.nz);
    end
end

jacConst_mp=[jacConst_mp;blz];
jacConst_mp=sparse(jacConst_mp);
end

 %------------- END OF CODE --------------                     






