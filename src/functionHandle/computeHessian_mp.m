function hessian_mp=computeHessian_mp(z_mp,sigma,lambda_mp,auxdata)
%computeHessian_mp - Evaluate the Hessian of the Lagrangian at z, for a
%multi-phase problem
%
% Syntax: hessian_mp=computeHessian_mp(z_mp,sigma,lambda_mp,auxdata)
%
% Inputs:
%    z - NLP variable 
%    sigma - Cost function factor
%    lambda - Vector of Lagrange multipliers
%    data - Structure with data needed to evaluate CostFn
%
% Outputs:
%    hessian - Sparse lower triangular segment of the hessian at z 
%
% Other m-files required: multipleShooting.m or directCollocation.m
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
hessian_mp=sparse(auxdata.mpdata.mpsizes.nz,auxdata.mpdata.mpsizes.nz);
nbl=auxdata.mpdata.mpsizes.nbl_nl;

x0=cell(auxdata.mpdata.mpsizes.nphase,1);xf=cell(auxdata.mpdata.mpsizes.nphase,1);
u0=cell(auxdata.mpdata.mpsizes.nphase,1);uf=cell(auxdata.mpdata.mpsizes.nphase,1);

for i=1:length(auxdata.phasedata)
    data=auxdata.phasedata{i};
    z=zeros(data.nz,1);
    z(data.zidx.org.z,1)=z_mp(data.zidx.mp.z,1);
    
    data.sigma=sigma;data.lambda=lambda_mp(data.zidx.mp.const,1);

    switch data.options.transcription

        case {'multiple_shooting'}

            hessian=multipleShooting('hessian',z,data);

        case {'direct_collocation'}

            switch data.options.discretization

                case{'globalLGR','hpLGR'} % p/hp Transcription Method

                    [hessian,XU0f]=directCollocationLGR('hessian',z,data,i);
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;
                     
                case {'resMinConstructionForSolution'}

                    hessian=resMinConstructionForSolution('hessian',z,data);   

                case {'resMinInterpolationForSolution'}

    %                 hessian=resMinInterpolationForSolution('hessian',z,data);   
                    [hessian,XU0f]=resMinOptimization('hessian',z,data);   
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;
                     
                otherwise % h Transcription Method

                    [hessian,XU0f]=directCollocation('hessian',z,data,i);
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;
                     
            end

        case {'integral_res_min'}

            [hessian,XU0f]=resMinOptimization('hessian',z,data);
             x0{i}=XU0f.x0;
             xf{i}=XU0f.xf;
             u0{i}=XU0f.u0;
             uf{i}=XU0f.uf;
                     
    end
    
    hessian_mp(data.hSidx.mp.XUXU.row,data.hSidx.mp.XUXU.col)=hessian_mp(data.hSidx.mp.XUXU.row,data.hSidx.mp.XUXU.col)+hessian(data.hSidx.org.XUXU.row,data.hSidx.org.XUXU.col);
    if size(hessian_mp(data.hSidx.org.PXU.row,data.hSidx.org.PXU.col),1) ~= size(hessian_mp(data.hSidx.mp.PXU.row,data.hSidx.mp.PXU.col),1)
        hessian_mp(data.hSidx.mp.PXU.row,data.hSidx.mp.PXU.col)=hessian_mp(data.hSidx.mp.PXU.row,data.hSidx.mp.PXU.col)+hessian(data.hSidx.org.PXU.row,data.hSidx.org.PXU.col)';
    else
        hessian_mp(data.hSidx.mp.PXU.row,data.hSidx.mp.PXU.col)=hessian_mp(data.hSidx.mp.PXU.row,data.hSidx.mp.PXU.col)+hessian(data.hSidx.org.PXU.row,data.hSidx.org.PXU.col);
    end
    hessian_mp(data.hSidx.mp.PP.row,data.hSidx.mp.PP.col)=hessian_mp(data.hSidx.mp.PP.row,data.hSidx.mp.PP.col)+hessian(data.hSidx.org.PP.row,data.hSidx.org.PP.col);
    if size(hessian_mp(data.hSidx.org.TXU.row,data.hSidx.org.TXU.col),1) ~= size(hessian_mp(data.hSidx.mp.TXU.row,data.hSidx.mp.TXU.col),1)
        hessian_mp(data.hSidx.mp.TXU.row,data.hSidx.mp.TXU.col)=hessian_mp(data.hSidx.mp.TXU.row,data.hSidx.mp.TXU.col)+hessian(data.hSidx.org.TXU.row,data.hSidx.org.TXU.col)';
    else
        hessian_mp(data.hSidx.mp.TXU.row,data.hSidx.mp.TXU.col)=hessian_mp(data.hSidx.mp.TXU.row,data.hSidx.mp.TXU.col)+hessian(data.hSidx.org.TXU.row,data.hSidx.org.TXU.col);
    end
    hessian_mp(data.hSidx.mp.TT.row,data.hSidx.mp.TT.col)=hessian_mp(data.hSidx.mp.TT.row,data.hSidx.mp.TT.col)+hessian(data.hSidx.org.TT.row,data.hSidx.org.TT.col);
%     
    
end


linkfunctions=auxdata.mpdata.linkfunctions;
mpdata=auxdata.mpdata;
t0=z_mp(end-mpdata.mpsizes.nt+1:end-1);
tf=z_mp(end-mpdata.mpsizes.nt+2:end);
p=z_mp(end-mpdata.mpsizes.nt-mpdata.mpsizes.np+1:end-mpdata.mpsizes.nt);
if isfield(data.data,'Pscale')
   p=scale_variables_back( p', data.data.Pscale, data.data.Pshift )';
end
    
bzz=spalloc(auxdata.mpdata.mpsizes.nz,auxdata.mpdata.mpsizes.nz,(auxdata.mpdata.mpsizes.nxu0f_inc(end)+auxdata.mpdata.mpsizes.nt+auxdata.mpdata.mpsizes.np).^2);
if nbl
idx=auxdata.mpdata.linkConst.all.idx;
nfd=size(idx,2);
mpdata=auxdata.mpdata;
adjoint=lambda_mp(end-nbl+1:end)';
e2=auxdata.mpdata.options.mp.perturbation.J.^2;

for i=1:nfd
   for j=1:i

    
    if j==i
        
         x0p=x0;x0m=x0;
         xfp=xf;xfm=xf;
         u0p=u0;u0m=u0;
         ufp=uf;ufm=uf;

        if mpdata.map.allidx(i)~=0
            x0p{mpdata.map.allidx(i)}=x0p{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.ex0(:,i);
            x0m{mpdata.map.allidx(j)}=x0m{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.ex0(:,j);
            xfp{mpdata.map.allidx(i)}=xfp{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.exf(:,i);
            xfm{mpdata.map.allidx(j)}=xfm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.exf(:,j);
            u0p{mpdata.map.allidx(i)}=u0p{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.eu0(:,i);
            u0m{mpdata.map.allidx(j)}=u0m{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.eu0(:,j);
            ufp{mpdata.map.allidx(i)}=ufp{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.euf(:,i);
            ufm{mpdata.map.allidx(j)}=ufm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.euf(:,j);
        end

         [~,blo]=linkfunctions(x0,xf,u0,uf,p,t0,tf,mpdata.data);
         [~,blp1]=linkfunctions(x0p,xfp,u0p,ufp,p+mpdata.map.ep(:,i),t0+mpdata.map.et0(:,i),tf+mpdata.map.etf(:,i),mpdata.data);
         [~,blp2]=linkfunctions(x0m,xfm,u0m,ufm,p-mpdata.map.ep(:,i),t0-mpdata.map.et0(:,i),tf-mpdata.map.etf(:,i),mpdata.data);
         blt=(blp2-2*blo+blp1).*adjoint'/e2; 
         bzz=bzz+sparse(idx(:,i),idx(:,j),blt,auxdata.mpdata.mpsizes.nz,auxdata.mpdata.mpsizes.nz); 
    else
         x0pp=x0;x0pm=x0;x0mp=x0;x0mm=x0;
         xfpp=xf;xfpm=xf;xfmp=xf;xfmm=xf;
         u0pp=u0;u0pm=u0;u0mp=u0;u0mm=u0;
         ufpp=uf;ufpm=uf;ufmp=uf;ufmm=uf;

        if mpdata.map.allidx(i)~=0
            x0pp{mpdata.map.allidx(i)}=x0pp{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.ex0(:,i);
            x0pp{mpdata.map.allidx(j)}=x0pp{mpdata.map.allidx(j)}+mpdata.map.phases{mpdata.map.allidx(j)}.ex0(:,j);

            x0pm{mpdata.map.allidx(i)}=x0pm{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.ex0(:,i);
            x0pm{mpdata.map.allidx(j)}=x0pm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.ex0(:,j);
            
            x0mp{mpdata.map.allidx(i)}=x0mp{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.ex0(:,i);
            x0mp{mpdata.map.allidx(j)}=x0mp{mpdata.map.allidx(j)}+mpdata.map.phases{mpdata.map.allidx(j)}.ex0(:,j);
            
            x0mm{mpdata.map.allidx(i)}=x0mm{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.ex0(:,i);
            x0mm{mpdata.map.allidx(j)}=x0mm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.ex0(:,j);
            
            xfpp{mpdata.map.allidx(i)}=xfpp{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.exf(:,i);
            xfpp{mpdata.map.allidx(j)}=xfpp{mpdata.map.allidx(j)}+mpdata.map.phases{mpdata.map.allidx(j)}.exf(:,j);

            xfpm{mpdata.map.allidx(i)}=xfpm{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.exf(:,i);
            xfpm{mpdata.map.allidx(j)}=xfpm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.exf(:,j);
            
            xfmp{mpdata.map.allidx(i)}=xfmp{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.exf(:,i);
            xfmp{mpdata.map.allidx(j)}=xfmp{mpdata.map.allidx(j)}+mpdata.map.phases{mpdata.map.allidx(j)}.exf(:,j);
            
            xfmm{mpdata.map.allidx(i)}=xfmm{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.exf(:,i);
            xfmm{mpdata.map.allidx(j)}=xfmm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.exf(:,j);
            
            u0pp{mpdata.map.allidx(i)}=u0pp{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.eu0(:,i);
            u0pp{mpdata.map.allidx(j)}=u0pp{mpdata.map.allidx(j)}+mpdata.map.phases{mpdata.map.allidx(j)}.eu0(:,j);

            u0pm{mpdata.map.allidx(i)}=u0pm{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.eu0(:,i);
            u0pm{mpdata.map.allidx(j)}=u0pm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.eu0(:,j);
            
            u0mp{mpdata.map.allidx(i)}=u0mp{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.eu0(:,i);
            u0mp{mpdata.map.allidx(j)}=u0mp{mpdata.map.allidx(j)}+mpdata.map.phases{mpdata.map.allidx(j)}.eu0(:,j);
            
            u0mm{mpdata.map.allidx(i)}=u0mm{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.eu0(:,i);
            u0mm{mpdata.map.allidx(j)}=u0mm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.eu0(:,j);
            
            ufpp{mpdata.map.allidx(i)}=ufpp{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.euf(:,i);
            ufpp{mpdata.map.allidx(j)}=ufpp{mpdata.map.allidx(j)}+mpdata.map.phases{mpdata.map.allidx(j)}.euf(:,j);

            ufpm{mpdata.map.allidx(i)}=ufpm{mpdata.map.allidx(i)}+mpdata.map.phases{mpdata.map.allidx(i)}.euf(:,i);
            ufpm{mpdata.map.allidx(j)}=ufpm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.euf(:,j);
            
            ufmp{mpdata.map.allidx(i)}=ufmp{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.euf(:,i);
            ufmp{mpdata.map.allidx(j)}=ufmp{mpdata.map.allidx(j)}+mpdata.map.phases{mpdata.map.allidx(j)}.euf(:,j);
            
            ufmm{mpdata.map.allidx(i)}=ufmm{mpdata.map.allidx(i)}-mpdata.map.phases{mpdata.map.allidx(i)}.euf(:,i);
            ufmm{mpdata.map.allidx(j)}=ufmm{mpdata.map.allidx(j)}-mpdata.map.phases{mpdata.map.allidx(j)}.euf(:,j);
        end
        
    [~,blpp]=linkfunctions(x0pp,xfpp,u0pp,ufpp,p+mpdata.map.ep(:,i)+mpdata.map.ep(:,j),t0+mpdata.map.et0(:,i)+mpdata.map.et0(:,j),tf+mpdata.map.etf(:,i)+mpdata.map.etf(:,j),mpdata.data);
    [~,blpm]=linkfunctions(x0pm,xfpm,u0pm,ufpm,p+mpdata.map.ep(:,i)-mpdata.map.ep(:,j),t0+mpdata.map.et0(:,i)-mpdata.map.et0(:,j),tf+mpdata.map.etf(:,i)-mpdata.map.etf(:,j),mpdata.data);
    [~,blmp]=linkfunctions(x0mp,xfmp,u0mp,ufmp,p-mpdata.map.ep(:,i)+mpdata.map.ep(:,j),t0-mpdata.map.et0(:,i)+mpdata.map.et0(:,j),tf-mpdata.map.etf(:,i)+mpdata.map.etf(:,j),mpdata.data);
    [~,blmm]=linkfunctions(x0mm,xfmm,u0mm,ufmm,p-mpdata.map.ep(:,i)-mpdata.map.ep(:,j),t0-mpdata.map.et0(:,i)-mpdata.map.et0(:,j),tf-mpdata.map.etf(:,i)-mpdata.map.etf(:,j),mpdata.data);
    blt=(blpp-blpm+blmm-blmp).*adjoint'/e2/4; 
   
    bzz=bzz+sparse(idx(:,i),idx(:,j),blt,auxdata.mpdata.mpsizes.nz,auxdata.mpdata.mpsizes.nz); 
    end
   end
end

end

hessian_mp=tril(sparse(hessian_mp)+bzz);
 %------------- END OF CODE --------------                     
