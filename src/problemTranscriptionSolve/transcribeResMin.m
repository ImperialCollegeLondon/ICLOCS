function [ auxdata ] = transcribeResMin(t_list,options,dataNLP )
%transcribeResErrorAnalysis - transcribe the problem in preparation for
%integrated residual minimization
%
% Syntax:  [ auxdata ] = transcribeResMin(t_list,options,dataNLP )
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

[nt,np,n,m,ng,~,M]=deal(dataNLP.sizes{1:7});


switch dataNLP.options.discretization
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        t0_list=t_list(1:end-1)'./t_list(end); tf_list=t_list(2:end)'./t_list(end); 
        if isfield(options,'pdegree') 
                npd=options.pdegree*ones(1,options.nsegment);
                npd_quad=max(max(npd)*4+1,5)*ones(1,length(t0_list));
                nps=options.nsegment;
        else
                npd=options.npsegment;
                npd_quad=max(max(npd)*4+1,5)*ones(1,length(t0_list));            
                nps=length(options.npsegment);                
        end
    case{'hermite'} 
        if M<5
            error('Direct interpolation with minimum residual not supported when the number of mesh nodes is less than 3 for Hermite-Simpson (h) method');
        end
        t0_list=t_list(1:end-1)'./t_list(end); tf_list=t_list(2:end)'./t_list(end); 
        npd=3*ones(1,length(t0_list));
        npd_quad=13*ones(1,length(t0_list));
        nps=length(t0_list);

    otherwise % h Transcription Method
        if strcmp(dataNLP.options.transcription,'integral_res_min')
            error('Direct interpolation with minimum residual not supported for the discretization method selected. Please use Hermite-Simpson (h) or LGR (p/hp) discretization method for direct collocation');
        else
            auxdata=[];
            return
        end
end

Mi=sum(npd);
npd=npd';npd_quad=npd_quad';
[npdu,~,npduidx]=unique([npd;npd_quad]);
idx_quad=find(npdu==npd_quad(1));
M_quad=sum(npd_quad+1); 

LGR=generateLGR_All(npdu); 

for i=1:length(npdu)
    LGR.interpWeights{1,i}=repmat([LGR.points{1,i};1],1,npdu(i)+1);
    LGR.interpWMat{1,i}=repmat(1./prod(LGR.interpWeights{1,i}-LGR.interpWeights{1,i}.'+eye(npdu(i)+1),1),npd_quad(1)+1,1);
    LGR.interpDist{1,i}=repmat([LGR.points{1,idx_quad};1],1,npdu(i)+1)-repmat([LGR.points{1,i};1].',npd_quad(1)+1,1);
    [LGR.interp_fixi{1,i},LGR.interp_fixj{1,i}]=find(LGR.interpDist{1,i}==0);
    LGR.interpDist{1,i}(LGR.interp_fixi{1,i},LGR.interp_fixj{1,i})=NaN;
    LGR.interpH{1,i}=LGR.interpWMat{1,i}./LGR.interpDist{1,i};
    LGR.sumInterpH{1,i}=sum(LGR.interpH{1,i},2);
end



        
switch dataNLP.options.discretization
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        DT=(tf_list-t0_list);
        if size(DT,1)>size(DT,2)
            DT=DT';
        end
        DT_seg_quad=repelem(DT',npd_quad+1,1);
        DT_seg_mesh=repelem(DT',npd,1);
        Hrowidxstart=1;
        Hcolidxstart=1;
        InterpH=zeros(M_quad,Mi);
        tau_quad=zeros(M_quad,1);
        sum_nps_quad=zeros(M_quad,1);
        sumInterpH=zeros(M_quad,1);
        interp_fixi=[];
        interp_fixj=[];
        cumsum_npd=[0;cumsum(npd)];
        cumsum_npd_quad=[0;cumsum(npd_quad+1)];
        for i=1:nps
              InterpH(Hrowidxstart:Hrowidxstart+npd_quad(i),Hcolidxstart:Hcolidxstart+npd(i))=LGR.interpH{1,npduidx(i)};
              interp_fixi=[interp_fixi;LGR.interp_fixi{1,npduidx(i)}+cumsum_npd_quad(i)];
              interp_fixj=[interp_fixj;LGR.interp_fixj{1,npduidx(i)}+cumsum_npd(i)];   
              tau_quad(Hrowidxstart:Hrowidxstart+npd_quad(i),1)=(tf_list(i)-t0_list(i))/2*[LGR.points{1,idx_quad};1]+(tf_list(i)+t0_list(i))/2;
              sum_nps_quad(Hrowidxstart:Hrowidxstart+npd_quad(i),i)=[LGR.weights{1, idx_quad};0];
              sumInterpH(Hrowidxstart:Hrowidxstart+npd_quad(i),1)=LGR.sumInterpH{1,npduidx(i)};
              Hrowidxstart=Hrowidxstart+npd_quad(i)+1;
              Hcolidxstart=Hcolidxstart+npd(i);  
        end
        D_mat=kron(speye(nps),sparse(LGR.diff_matrix_square{idx_quad}));

        tau_inc=dataNLP.tau_inc;
        tf=tf_list(end);
        t0=t0_list(1);
        tau=(tf-t0)/2*tau_inc+(tf+t0)/2;
        
        idx_perturb=zeros(length(dataNLP.tau_order_idx)+1,max(dataNLP.tau_order_idx)+1);
        dataNLP.tau_order_idx=[dataNLP.tau_order_idx;1];
        for i=1:max(dataNLP.tau_order_idx)
            if i==1
                idx_mesh=find(dataNLP.tau_order_idx==1);
                idx_mesh_1=idx_mesh(1:2:end);
                idx_mesh_2=idx_mesh(2:2:end);
                
                pt_mesh=dataNLP.tau_order_idx==1;
                pt_mesh_1=pt_mesh;pt_mesh_2=pt_mesh;
                pt_mesh_1(idx_mesh_2)=0;
                pt_mesh_2(idx_mesh_1)=0;
                
                idx_perturb(:,1)=pt_mesh_1;
                idx_perturb(:,2)=pt_mesh_2;
            else
                idx_perturb(:,i+1)=dataNLP.tau_order_idx==i;
            end
        end
        
        idx_perturb_hes=zeros(length(dataNLP.tau_order_idx_hes),max(npd)*3);
        idx_perturb_hes(:,1)=dataNLP.tau_order_idx_hes==1;
        idx_perturb_hes(:,2)=dataNLP.tau_order_idx_hes==101;
        idx_perturb_hes(:,3)=dataNLP.tau_order_idx_hes==201;
        for i=2:max(npd)
            idx_perturb_hes(:,2+i)=dataNLP.tau_order_idx_hes==i;
            idx_perturb_hes(:,2+max(npd)-1+i)=dataNLP.tau_order_idx_hes==i+100;
            idx_perturb_hes(:,2+(max(npd)-1)*2+i)=dataNLP.tau_order_idx_hes==i+200;
        end

        InterpH(isnan(InterpH))=0;
        sumInterpH(isnan(sumInterpH))=0;
        
        InterpH=sparse(InterpH);
        auxdata.InterpH=InterpH;
        auxdata.sumInterpH=sumInterpH;
        auxdata.sumInterpHMat=sparse(diag(1./sumInterpH));
        auxdata.interp_fixi=interp_fixi;
        auxdata.interp_fixj=interp_fixj;
        auxdata.D_mat=sparse(D_mat);
        auxdata.DT_seg=DT_seg_quad;
        auxdata.DT_seg_mat_d2=diag(1./(auxdata.DT_seg/2));
        auxdata.sum_nps_quad=sparse(sum_nps_quad');
        auxdata.idx_perturb=idx_perturb;
        auxdata.idx_perturb_hes=idx_perturb_hes;
        
        auxdata.free_time=1;
        
        if dataNLP.options.scaling
            auxdata.ResNormScale=dataNLP.data.Xscale.*dataNLP.data.resNormCusWeight;
            ResNormScaleMat=repmat(auxdata.ResNormScale , nps, 1 );
            auxdata.ResNormScaleMat=sparse(diag(ResNormScaleMat(:)));
            if isfield(dataNLP.data,'discErrorConstScaling')
                ResConstScaleMat=repmat(dataNLP.data.discErrorConstScaling , nps, 1 );
                auxdata.ResConstScaleMat=sparse(diag(ResConstScaleMat(:)));
            end
        else
            auxdata.ResNormScale=dataNLP.data.resNormCusWeight;
            ResNormScaleMat=repmat( auxdata.ResNormScale, nps, 1 );
            auxdata.ResNormScaleMat=sparse(diag(ResNormScaleMat(:)));
            if isfield(dataNLP.data,'discErrorConstScaling')
                ResConstScaleMat=repmat(dataNLP.data.discErrorConstScaling , nps, 1 );
                auxdata.ResConstScaleMat=sparse(diag(ResConstScaleMat(:)));
            end
        end

        
    case{'hermite'} 
        DT=(tf_list-t0_list);
        DT_seg_quad=repelem(DT',npd_quad+1,1);
        DxHS1=sparse(diag(-3*ones(M,1))+diag(4*ones(M-1,1),1)+diag(-ones(M-2,1),2));
        DxHS1([2:2:end,end],:)=[];
        DxHS2=sparse(diag(-ones(M,1))+diag(ones(M-2,1),2));
        DxHS2([2:2:end,end],:)=[];
        DxHS3=sparse(diag(1*ones(M,1))+diag(-4*ones(M-1,1),1)+diag(3*ones(M-2,1),2));
        DxHS3([2:2:end,end],:)=[];
        
        DxHS_hf=sparse(diag(-5/2*ones(M,1))+diag(2*ones(M-1,1),1)+diag(0.5*ones(M-2,1),2));
        DxHS_hf([2:2:end,end],:)=[];
        DxHS_hf=DxHS_hf./DT';
        DxHS_p1=sparse(diag(4*ones(M,1))+diag(-8*ones(M-1,1),1)+diag(4*ones(M-2,1),2));
        DxHS_p1([2:2:end,end],:)=[];
        DxHS_p1=DxHS_p1./DT';
        
        DT_seg_mesh=repelem(DT',npd,1);
        
        tau_quad=DT_seg_quad/2.*repmat([LGR.points{1,idx_quad};1],nps,1)+repelem(tf_list'+t0_list',npd_quad+1,1)./2;
        DT_tau=tau_quad-repelem(t0_list',npd_quad+1,1);
        u_coef1=2./((DT_seg_quad).^2).*(DT_tau-DT_seg_quad/2).*(DT_tau-DT_seg_quad);
        u_coef2=-4./((DT_seg_quad).^2).*DT_tau.*(DT_tau-DT_seg_quad);
        u_coef3=2./((DT_seg_quad).^2).*(DT_tau-DT_seg_quad/2).*DT_tau;
        
        x_coef1=DT_tau+0.5*(-3)*DT_tau.^2./DT_seg_quad+1/3*2*DT_tau.^3./DT_seg_quad.^2;
        x_coef2=0.5*4*DT_tau.^2./DT_seg_quad+1/3*(-4)*DT_tau.^3./DT_seg_quad.^2;
        x_coef3=-0.5*DT_tau.^2./DT_seg_quad+1/3*2*DT_tau.^3./DT_seg_quad.^2;
        
        AuHS=zeros(M_quad,M);
        AfHS=zeros(M_quad,M);
        sum_nps_quad=zeros(M_quad,1);
        AxHS=zeros(M_quad,3*(nps));
        AHSrowidxstart=1;
        AHScolidxstart=1;
        AuHScolidxstart=1;
        for i=1:nps
            AfHS(AHSrowidxstart:AHSrowidxstart+npd_quad(i),AHScolidxstart:AHScolidxstart+2)=[u_coef1(AHSrowidxstart:AHSrowidxstart++npd_quad(i)) u_coef2(AHSrowidxstart:AHSrowidxstart++npd_quad(i)) u_coef3(AHSrowidxstart:AHSrowidxstart++npd_quad(i))];
            AuHS(AHSrowidxstart:AHSrowidxstart+npd_quad(i),AuHScolidxstart:AuHScolidxstart+2)=[u_coef1(AHSrowidxstart:AHSrowidxstart++npd_quad(i)) u_coef2(AHSrowidxstart:AHSrowidxstart++npd_quad(i)) u_coef3(AHSrowidxstart:AHSrowidxstart++npd_quad(i))];
            AxHS(AHSrowidxstart:AHSrowidxstart+npd_quad(i),AHScolidxstart:AHScolidxstart+2)=[x_coef1(AHSrowidxstart:AHSrowidxstart++npd_quad(i)) x_coef2(AHSrowidxstart:AHSrowidxstart++npd_quad(i)) x_coef3(AHSrowidxstart:AHSrowidxstart++npd_quad(i))];
            sum_nps_quad(AHSrowidxstart:AHSrowidxstart+npd_quad(i),i)=[LGR.weights{1, idx_quad};0];
            AHSrowidxstart=AHSrowidxstart+npd_quad(i)+1;
            AHScolidxstart=AHScolidxstart+3;
            AuHScolidxstart=AuHScolidxstart+2;
        end
        
        repXend_mat=eye((M-1)/2);
        repXend_mat=repelem(repXend_mat,npd_quad+1,1);
        

        AuHS=sparse(AuHS);
        AfHS=sparse(AfHS);
        AxHS=sparse(AxHS);
        
        tau=[0;dataNLP.tau_inc]/2;
        
        idx_perturb=zeros(M,3);
        idx_perturb(1:4:end,1)=1;
        idx_perturb(3:4:end,2)=1;
        idx_perturb(2:2:end,3)=1;
        
        idx_perturb_hes=zeros(M,6);
        idx_perturb_hes(1:6:end,1)=1;
        idx_perturb_hes(3:6:end,2)=1;
        idx_perturb_hes(5:6:end,3)=1;
        idx_perturb_hes(2:6:end,4)=1;
        idx_perturb_hes(4:6:end,5)=1;
        idx_perturb_hes(6:6:end,6)=1;
        
        auxdata.DxHS1=sparse(DxHS1);
        auxdata.DxHS2=sparse(DxHS2);
        auxdata.DxHS3=sparse(DxHS3);
        auxdata.DxHS_hf=sparse(DxHS_hf);
        auxdata.DxHS_p1=sparse(DxHS_p1);
        
        auxdata.idxcoll=eye(M);
        auxdata.idxcoll(2:2:end,:)=[];
        auxdata.idx_m1=eye((M+1)/2);
        auxdata.idx_m1(end,:)=[];

        auxdata.AuHS=AuHS;
        auxdata.AxHS=AxHS;
        auxdata.AfHS=AfHS;
        auxdata.repXend_mat=sparse(repXend_mat);
        auxdata.sum_nps_quad=sparse(sum_nps_quad');
        auxdata.idx_perturb=idx_perturb;
        auxdata.idx_perturb_hes=idx_perturb_hes;
        
        if dataNLP.options.scaling

            auxdata.ResNormScale=dataNLP.data.Xscale.*dataNLP.data.resNormCusWeight;
            ResNormScaleMat=repmat(auxdata.ResNormScale , nps, 1 );
            auxdata.ResNormScaleMat=sparse(diag(ResNormScaleMat(:)));
            if isfield(dataNLP.data,'discErrorConstScaling')
                ResConstScaleMat=repmat(dataNLP.data.discErrorConstScaling , nps, 1 );
                auxdata.ResConstScaleMat=sparse(diag(ResConstScaleMat(:)));
            end
        else
            auxdata.ResNormScale=dataNLP.data.resNormCusWeight;
            ResNormScaleMat=repmat( auxdata.ResNormScale, nps, 1 );
            auxdata.ResNormScaleMat=sparse(diag(ResNormScaleMat(:)));
            if isfield(dataNLP.data,'discErrorConstScaling')
                ResConstScaleMat=repmat(dataNLP.data.discErrorConstScaling , nps, 1 );
                auxdata.ResConstScaleMat=sparse(diag(ResConstScaleMat(:)));
            end
        end

        
        
end

idx_spt=1:(npd_quad(1)+1):M_quad;
idx_ept=(npd_quad(1)+1):(npd_quad(1)+1):M_quad;

auxdata.t0_list=t0_list;auxdata.tf_list=tf_list;
auxdata.LGR=LGR;auxdata.idx_quad=idx_quad;
auxdata.nps=nps;auxdata.npd=npd;

auxdata.tau_quad=tau_quad;
auxdata.M_quad=M_quad;
auxdata.DT_seg_quad=DT_seg_quad;
auxdata.DT_seg_mesh=DT_seg_mesh;
auxdata.tau=tau;
        
auxdata.idx_spt=idx_spt;
auxdata.idx_ept=idx_ept;

auxdata.dataNLP=dataNLP;
auxdata.npd_quad=npd_quad;


end

