function [infoNLP,auxdata,options]=transcribeMultiphaseOCP(problem,guess,options)
%transcribeMultiphaseOCP - Process information from 'problem', 'guess' and
%'options' for NLP solver, for multiphase problems
% 
%Specifically:
%Error checking of function definitions and bounds
%Define bounds for NLP variable + continuity, path and boundary constraints
%Format matrices for direct transcription method(if required)
%Generate initial guess for optimization
%Generate structure of the jacobian of the constraints(if required)
%Construct optimal finite-difference pertubation sets(if required)
%
% Syntax:  [infoNLP,auxdata,options]=transcribeMultiphaseOCP(problem,guess,options)
%
% Inputs:
%    problem - Optimal control problem definition
%    guess   - Guess used to generate starting point for optimization
%    options - Settings in file settings.m
%
% Outputs:
%    infoNLP - Information required by the NLP solver
%    auxdata - Data passed to the functions evaluated during optimization
%    options - Settings after processing
%
% Other m-files required: scale_problem, getStructure, getStructureA, getPertubations.
% Subfunctions: checkErrors, transcriptionMatrix
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


global sol;  
sol=[];                             % Initialize solution structure

if strcmp(options.mp.transcription,'integral_res_min')
    error('Multi-phase problem currently only support Direct Collcation transcription')
end
    
% Running mode
if ~isfield(problem.mp.data,'mode')
    problem.mp.data.mode.currentMode='Original';
end

if isempty(options.mp.perturbation.J)
  options.mp.perturbation.J=(eps/2)^(1/3);
end
if isempty(options.mp.perturbation.H)    
  options.mp.perturbation.H=(8*eps)^(1/3);
end    



problem.mp.data.transcription=options.mp.transcription;
data.data=problem.mp.data;

infoNLP_XU.zl=[];
infoNLP_XU.zu=[];
infoNLP_XU.z0=[];

infoNLP_XUg.cl=[];
infoNLP_XUg.cu=[];
infoNLP_b.cl=[];
infoNLP_b.cu=[];


phaseinfoNLP=cell(1,length(problem.phases));
phasedata=cell(1,length(problem.phases));
data.linkConst.xu0f.idx=[];
if options.mp.scaling
    data.linkConst.xu0f.scaling.scale=[];
    data.linkConst.xu0f.scaling.shift=[];
end
data.mpsizes.nzblock=[];
data.mpsizes.nzblock_inc=[];
nzblock_current=0;
data.mpsizes.nxu0f=[];

data.mpsizes.nphase=length(problem.phases);
for i=1:data.mpsizes.nphase
    
    
    
    
    if ~isfield(options.phaseoptions{i},'transcription')
        names = [fieldnames(options.mp); fieldnames(options.phaseoptions{i})];
        options.phaseoptions{i} = cell2struct([struct2cell(options.mp); struct2cell(options.phaseoptions{i})], names, 1);
        
    end
    
    options.phaseoptions{i}.perturbation.J=options.mp.perturbation.J;
    options.phaseoptions{i}.perturbation.H=options.mp.perturbation.H;
    [phaseinfoNLP{i},phasedata{i},options.phaseoptions{i}]=transcribeOCP_eachPhase(problem.phases{i},guess.phases{i},options.phaseoptions{i});
    

    phasedata{i}.data.mode.currentMode=problem.mp.data.mode.currentMode;
    
    if strcmp(options.mp.transcription,'integral_res_min')
        data_phase=phasedata{i}.dataNLP;
    else
        data_phase=phasedata{i};
    end
    
    data.linkConst.xu0f.idx=[data.linkConst.xu0f.idx [data_phase.infoForLinkConst.x0idx data_phase.infoForLinkConst.xfidx data_phase.infoForLinkConst.u0idx data_phase.infoForLinkConst.ufidx]+nzblock_current];
    if options.mp.scaling
        data.linkConst.xu0f.scaling.scale=[data.linkConst.xu0f.scaling.scale [data_phase.data.Xscale data_phase.data.Xscale data_phase.data.Uscale data_phase.data.Uscale]];
        data.linkConst.xu0f.scaling.shift=[data.linkConst.xu0f.scaling.shift [data_phase.data.Xshift data_phase.data.Xshift data_phase.data.Ushift data_phase.data.Ushift]];
    end
    data_phase.hSidx.mp.XUXU.row=length(infoNLP_XU.zl)+1:length(infoNLP_XU.zl)+length(data_phase.hSidx.org.XUXU.row);
    data_phase.hSidx.mp.XUXU.col=data_phase.hSidx.mp.XUXU.row;
    data_phase.zidx.mp.xu=data_phase.hSidx.mp.XUXU.row;
    infoNLP_XU.zl=[infoNLP_XU.zl;phaseinfoNLP{i}.zl(data_phase.hSidx.org.XUXU.row,1)];
    infoNLP_XU.zu=[infoNLP_XU.zu;phaseinfoNLP{i}.zu(data_phase.hSidx.org.XUXU.row,1)];
    infoNLP_XU.z0=[infoNLP_XU.z0;phaseinfoNLP{i}.z0(data_phase.hSidx.org.XUXU.row,1)];
    
    nzblock_current=length(infoNLP_XU.zl);
    data.mpsizes.nzblock_inc=[data.mpsizes.nzblock_inc length(infoNLP_XU.zl)];
    data.mpsizes.nzblock=[data.mpsizes.nzblock length(data_phase.hSidx.org.XUXU.row)];
    data.mpsizes.nxu0f=[data.mpsizes.nxu0f data_phase.infoForLinkConst.nxu0f];
    
    data_phase.jSidx.mp.XUg.row=length(infoNLP_XUg.cl)+1:length(infoNLP_XUg.cl)+length(data_phase.jSidx.org.XUg.row);
    data_phase.jSidx.mp.XUg.col=data_phase.hSidx.mp.XUXU.col;
    
    
    infoNLP_XUg.cl=[infoNLP_XUg.cl;phaseinfoNLP{i}.cl(data_phase.jSidx.org.XUg.row,1)];
    infoNLP_XUg.cu=[infoNLP_XUg.cu;phaseinfoNLP{i}.cu(data_phase.jSidx.org.XUg.row,1)];
    infoNLP_b.cl=[infoNLP_b.cl;phaseinfoNLP{i}.cl(data_phase.jSidx.org.B_XUg.row,1)];
    infoNLP_b.cu=[infoNLP_b.cu;phaseinfoNLP{i}.cu(data_phase.jSidx.org.B_XUg.row,1)];
    
    if strcmp(options.mp.transcription,'integral_res_min')
        phasedata{i}.dataNLP=data_phase;
    else
        phasedata{i}=data_phase;
    end
end 

data.mpsizes.nxu=length(infoNLP_XU.zl);
data.mpsizes.np=length(problem.mp.parameters.pl);
data.mpsizes.nt=length(problem.mp.time.t_min);
data.mpsizes.nxug=length(infoNLP_XUg.cl);
data.mpsizes.nb=length(infoNLP_b.cl);
data.mpsizes.nbl_l=length(problem.mp.constraints.bll.linear);
data.mpsizes.nbl_nl=length(problem.mp.constraints.bll.nonlinear);
data.mpsizes.nbl=data.mpsizes.nbl_l+data.mpsizes.nbl_nl;
data.mpsizes.nConst=data.mpsizes.nxug+data.mpsizes.nb+data.mpsizes.nbl;
data.mpsizes.nz=data.mpsizes.nxu+data.mpsizes.np+data.mpsizes.nt;
data.mpsizes.nxu0f_inc=cumsum([0 data.mpsizes.nxu0f]);
data.mpsizes.nxupt0f=sum(data.mpsizes.nxu0f)+data.mpsizes.np+data.mpsizes.nt;
data.linkConst.p.idx=data.mpsizes.nzblock_inc(end)+1:data.mpsizes.nzblock_inc(end)+data.mpsizes.np;
data.linkConst.t.idx=data.mpsizes.nzblock_inc(end)+data.mpsizes.np+1:data.mpsizes.nzblock_inc(end)+data.mpsizes.np+data.mpsizes.nt;

map=speye(data.mpsizes.nxupt0f);
map_ez=options.mp.perturbation.J.*map;
map_idx=zeros(1,data.mpsizes.nxupt0f);
data.map.all=map;
data.map.ez=map_ez;
data.map.phases=cell(length(problem.phases),1);
data.map.et0=map_ez((data.mpsizes.nxupt0f-data.mpsizes.nt+1):data.mpsizes.nxupt0f-1,:);
data.map.etf=map_ez((data.mpsizes.nxupt0f-data.mpsizes.nt+2):data.mpsizes.nxupt0f,:);
if options.mp.scaling
    if strcmp(options.mp.transcription,'integral_res_min')
        if isfield(phasedata{1}.dataNLP.data,'Pscale')
            data.map.ep=map_ez((data.mpsizes.nxupt0f-data.mpsizes.nt-data.mpsizes.np+1):(data.mpsizes.nxupt0f-data.mpsizes.nt),:)./phasedata{1}.dataNLP.data.Pscale';
            infoNLP_P.zl=scale_variables( problem.mp.parameters.pl, phasedata{1}.dataNLP.data.Pscale, phasedata{1}.dataNLP.data.Pshift )';
            infoNLP_P.zu=scale_variables( problem.mp.parameters.pu, phasedata{1}.dataNLP.data.Pscale, phasedata{1}.dataNLP.data.Pshift )';
            infoNLP_P.z0=scale_variables( guess.mp.parameters, phasedata{1}.dataNLP.data.Pscale, phasedata{1}.dataNLP.data.Pshift )';
        else
            data.map.ep=map_ez((data.mpsizes.nxupt0f-data.mpsizes.nt-data.mpsizes.np+1):(data.mpsizes.nxupt0f-data.mpsizes.nt),:);
            infoNLP_P.zl=[];
            infoNLP_P.zu=[];
            infoNLP_P.z0=[];
        end
    else
        if isfield(phasedata{1}.data,'Pscale')
            data.map.ep=map_ez((data.mpsizes.nxupt0f-data.mpsizes.nt-data.mpsizes.np+1):(data.mpsizes.nxupt0f-data.mpsizes.nt),:)./phasedata{1}.data.Pscale';
            infoNLP_P.zl=scale_variables( problem.mp.parameters.pl, phasedata{1}.data.Pscale, phasedata{1}.data.Pshift )';
            infoNLP_P.zu=scale_variables( problem.mp.parameters.pu, phasedata{1}.data.Pscale, phasedata{1}.data.Pshift )';
            infoNLP_P.z0=scale_variables( guess.mp.parameters, phasedata{1}.data.Pscale, phasedata{1}.data.Pshift )';
        else
            data.map.ep=map_ez((data.mpsizes.nxupt0f-data.mpsizes.nt-data.mpsizes.np+1):(data.mpsizes.nxupt0f-data.mpsizes.nt),:);
            infoNLP_P.zl=[];
            infoNLP_P.zu=[];
            infoNLP_P.z0=[];
        end
    end
else
    data.map.ep=map_ez((data.mpsizes.nxupt0f-data.mpsizes.nt-data.mpsizes.np+1):(data.mpsizes.nxupt0f-data.mpsizes.nt),:);
    infoNLP_P.zl=problem.mp.parameters.pl';
    infoNLP_P.zu=problem.mp.parameters.pu';
    infoNLP_P.z0=guess.mp.parameters';
end


nxupt0f=cumsum([0 data.mpsizes.nxu0f]);
nbmp=0;
for i=1:data.mpsizes.nphase
    
    if strcmp(options.mp.transcription,'integral_res_min')
        data_phase=phasedata{i}.dataNLP;
    else
        data_phase=phasedata{i};
    end
    
    phase_n=data_phase.sizes{3};
    phase_m=data_phase.sizes{4};
    if options.mp.scaling
        data.map.phases{i}.ex0=map_ez(nxupt0f(i)+1:nxupt0f(i)+phase_n,:)./data_phase.data.Xscale';
        data.map.phases{i}.exf=map_ez(nxupt0f(i)+phase_n+1:nxupt0f(i)+phase_n*2,:)./data_phase.data.Xscale';
        data.map.phases{i}.eu0=map_ez(nxupt0f(i)+phase_n*2+1:nxupt0f(i)+phase_n*2+phase_m,:)./data_phase.data.Uscale';
        data.map.phases{i}.euf=map_ez(nxupt0f(i)+phase_n*2+phase_m+1:nxupt0f(i)+phase_n*2+phase_m*2,:)./data_phase.data.Uscale';
    else
        data.map.phases{i}.ex0=map_ez(nxupt0f(i)+1:nxupt0f(i)+phase_n,:);
        data.map.phases{i}.exf=map_ez(nxupt0f(i)+phase_n+1:nxupt0f(i)+phase_n*2,:);
        data.map.phases{i}.eu0=map_ez(nxupt0f(i)+phase_n*2+1:nxupt0f(i)+phase_n*2+phase_m,:);
        data.map.phases{i}.euf=map_ez(nxupt0f(i)+phase_n*2+phase_m+1:nxupt0f(i)+phase_n*2+phase_m*2,:);
    end
    map_idx(nxupt0f(i)+1:nxupt0f(i)+phase_n*2+phase_m*2)=i;
    
    if data.mpsizes.np
        data_phase.hSidx.mp.PXU.row=data.mpsizes.nxu+1:data.mpsizes.nxu+data.mpsizes.np;
        data_phase.hSidx.mp.PXU.col=data_phase.hSidx.mp.XUXU.col;
        data_phase.hSidx.mp.PP.row=data.mpsizes.nxu+1:data.mpsizes.nxu+data.mpsizes.np;
        data_phase.hSidx.mp.PP.col=data_phase.hSidx.mp.PP.row;
        data_phase.zidx.mp.p=data_phase.hSidx.mp.PP.row;
        
        data_phase.jSidx.mp.P.row=data_phase.jSidx.mp.XUg.row;
        data_phase.jSidx.mp.P.col=data.mpsizes.nxu+1:data.mpsizes.nxu+data.mpsizes.np;
    else
        data_phase.hSidx.mp.PXU.row=[];
        data_phase.hSidx.mp.PXU.col=[];
        data_phase.hSidx.mp.PP.row=[];
        data_phase.hSidx.mp.PP.col=[];
        data_phase.zidx.mp.p=[];
        
        data_phase.jSidx.mp.P.row=[];
        data_phase.jSidx.mp.P.col=[];
    end
    
    nt=data_phase.sizes{1};
    nb=data_phase.sizes{6};
    if nt
        data_phase.hSidx.mp.TXU.row=[data.mpsizes.nxu+data.mpsizes.np+problem.phases{i}.time.t0_idx, data.mpsizes.nxu+data.mpsizes.np+problem.phases{i}.time.tf_idx];
        data_phase.hSidx.mp.TXU.col=data_phase.hSidx.mp.XUXU.col;
        data_phase.hSidx.mp.TT.row=[data.mpsizes.nxu+data.mpsizes.np+problem.phases{i}.time.t0_idx, data.mpsizes.nxu+data.mpsizes.np+problem.phases{i}.time.tf_idx];
        data_phase.hSidx.mp.TT.col=data_phase.hSidx.mp.TT.row;
        data_phase.zidx.mp.t=data_phase.hSidx.mp.TT.row;
        
        data_phase.jSidx.mp.T.row=data_phase.jSidx.mp.XUg.row;
        data_phase.jSidx.mp.T.col=[data.mpsizes.nxu+data.mpsizes.np+problem.phases{i}.time.t0_idx, data.mpsizes.nxu+data.mpsizes.np+problem.phases{i}.time.tf_idx];
    else
        data_phase.hSidx.mp.TXU.row=[];
        data_phase.hSidx.mp.TXU.col=[];
        data_phase.hSidx.mp.TT.row=[];
        data_phase.hSidx.mp.TT.col=[];
        data_phase.zidx.mp.t=[];
        
        data_phase.jSidx.mp.T.row=[];
        data_phase.jSidx.mp.T.col=[];
    end
        
    if nb
        data_phase.jSidx.mp.B_XUg.row=data.mpsizes.nxug+nbmp+1:data.mpsizes.nxug+nbmp+length(data_phase.jSidx.org.B_XUg.row);
        data_phase.jSidx.mp.B_XUg.col=data_phase.hSidx.mp.XUXU.col;
        if data.mpsizes.np
            data_phase.jSidx.mp.B_P.row=data.mpsizes.nxug+nbmp+1:data.mpsizes.nxug+nbmp+length(data_phase.jSidx.org.B_XUg.row);
            data_phase.jSidx.mp.B_P.col=data_phase.jSidx.mp.P.col;
        else
            data_phase.jSidx.mp.B_P.row=[];
            data_phase.jSidx.mp.B_P.col=[];
        end
        if nt
            data_phase.jSidx.mp.B_T.row=data.mpsizes.nxug+nbmp+1:data.mpsizes.nxug+nbmp+length(data_phase.jSidx.org.B_XUg.row);
            data_phase.jSidx.mp.B_T.col=data_phase.jSidx.mp.T.col;
        else
            data_phase.jSidx.mp.B_T.row=[];
            data_phase.jSidx.mp.B_T.col=[];
        end
        nbmp=nbmp+length(data_phase.jSidx.org.B_XUg.row);
    else
        data_phase.jSidx.mp.B_XUg.row=[];
        data_phase.jSidx.mp.B_XUg.col=[];
        data_phase.jSidx.mp.B_P.row=[];
        data_phase.jSidx.mp.B_P.col=[];
        data_phase.jSidx.mp.B_T.row=[];
        data_phase.jSidx.mp.B_T.col=[];
    end
    
    data_phase.zidx.mp.z=[data_phase.zidx.mp.t data_phase.zidx.mp.p data_phase.zidx.mp.xu];
    data_phase.zidx.org.z=[data_phase.zidx.org.t data_phase.zidx.org.p data_phase.zidx.org.xu];
    
    data_phase.zidx.mp.const=[data_phase.jSidx.mp.XUg.row data_phase.jSidx.mp.B_XUg.row];
    data_phase.zidx.org.const=[data_phase.jSidx.org.XUg.row data_phase.jSidx.org.B_XUg.row];
    
    data_phase.jSidx.mp.all.row=[data_phase.jSidx.mp.XUg.row data_phase.jSidx.mp.P.row data_phase.jSidx.mp.T.row data_phase.jSidx.mp.B_XUg.row data_phase.jSidx.mp.B_P.row data_phase.jSidx.mp.B_T.row];
    data_phase.jSidx.mp.all.col=[data_phase.jSidx.mp.XUg.col data_phase.jSidx.mp.P.col data_phase.jSidx.mp.T.col data_phase.jSidx.mp.B_XUg.col data_phase.jSidx.mp.B_P.col data_phase.jSidx.mp.B_T.col];
    data_phase.jSidx.org.all.row=[data_phase.jSidx.org.XUg.row data_phase.jSidx.org.P.row data_phase.jSidx.org.T.row data_phase.jSidx.org.B_XUg.row data_phase.jSidx.org.B_P.row data_phase.jSidx.org.B_T.row];
    data_phase.jSidx.org.all.col=[data_phase.jSidx.org.XUg.col data_phase.jSidx.org.P.col data_phase.jSidx.org.T.col data_phase.jSidx.org.B_XUg.col data_phase.jSidx.org.B_P.col data_phase.jSidx.org.B_T.col];

    [data_phase.jSidx.mp.all.rowidx,data_phase.jSidx.mp.all.colidx]=meshgrid(data_phase.jSidx.mp.all.row,data_phase.jSidx.mp.all.col);
    
    if strcmp(options.mp.transcription,'integral_res_min')
        phasedata{i}.dataNLP=data_phase;
    else
        phasedata{i}=data_phase;
    end
end 
data.map.allidx=map_idx;

data.hessianStruct=sparse(data.mpsizes.nz,data.mpsizes.nz);
data.jacStruct=sparse(data.mpsizes.nConst,data.mpsizes.nz);
for i=1:data.mpsizes.nphase
    
    if strcmp(options.mp.transcription,'integral_res_min')
        data_phase=phasedata{i}.dataNLP;
    else
        data_phase=phasedata{i};
    end
    
    data.hessianStruct(data_phase.hSidx.mp.XUXU.row,data_phase.hSidx.mp.XUXU.col)=data.hessianStruct(data_phase.hSidx.mp.XUXU.row,data_phase.hSidx.mp.XUXU.col)+data_phase.hessianStruct(data_phase.hSidx.org.XUXU.row,data_phase.hSidx.org.XUXU.col);
    if size(data_phase.hessianStruct(data_phase.hSidx.org.PXU.row,data_phase.hSidx.org.PXU.col),1) ~= size(data.hessianStruct(data_phase.hSidx.mp.PXU.row,data_phase.hSidx.mp.PXU.col),1)
        data.hessianStruct(data_phase.hSidx.mp.PXU.row,data_phase.hSidx.mp.PXU.col)=data.hessianStruct(data_phase.hSidx.mp.PXU.row,data_phase.hSidx.mp.PXU.col)+data_phase.hessianStruct(data_phase.hSidx.org.PXU.row,data_phase.hSidx.org.PXU.col)';
    else
        data.hessianStruct(data_phase.hSidx.mp.PXU.row,data_phase.hSidx.mp.PXU.col)=data.hessianStruct(data_phase.hSidx.mp.PXU.row,data_phase.hSidx.mp.PXU.col)+data_phase.hessianStruct(data_phase.hSidx.org.PXU.row,data_phase.hSidx.org.PXU.col);
    end
    data.hessianStruct(data_phase.hSidx.mp.PP.row,data_phase.hSidx.mp.PP.col)=data.hessianStruct(data_phase.hSidx.mp.PP.row,data_phase.hSidx.mp.PP.col)+data_phase.hessianStruct(data_phase.hSidx.org.PP.row,data_phase.hSidx.org.PP.col);
    if size(data_phase.hessianStruct(data_phase.hSidx.org.TXU.row,data_phase.hSidx.org.TXU.col),1) ~= size(data.hessianStruct(data_phase.hSidx.mp.TXU.row,data_phase.hSidx.mp.TXU.col),1)
        data.hessianStruct(data_phase.hSidx.mp.TXU.row,data_phase.hSidx.mp.TXU.col)=data.hessianStruct(data_phase.hSidx.mp.TXU.row,data_phase.hSidx.mp.TXU.col)+data_phase.hessianStruct(data_phase.hSidx.org.TXU.row,data_phase.hSidx.org.TXU.col)';
    else
        data.hessianStruct(data_phase.hSidx.mp.TXU.row,data_phase.hSidx.mp.TXU.col)=data.hessianStruct(data_phase.hSidx.mp.TXU.row,data_phase.hSidx.mp.TXU.col)+data_phase.hessianStruct(data_phase.hSidx.org.TXU.row,data_phase.hSidx.org.TXU.col);
    end
    data.hessianStruct(data_phase.hSidx.mp.TT.row,data_phase.hSidx.mp.TT.col)=data.hessianStruct(data_phase.hSidx.mp.TT.row,data_phase.hSidx.mp.TT.col)+data_phase.hessianStruct(data_phase.hSidx.org.TT.row,data_phase.hSidx.org.TT.col);
    data.jacStruct(data_phase.jSidx.mp.all.row,data_phase.jSidx.mp.all.col)=data.jacStruct(data_phase.jSidx.mp.all.row,data_phase.jSidx.mp.all.col)+data_phase.jacStruct(data_phase.jSidx.org.all.row,data_phase.jSidx.org.all.col);

end
data.jacStruct(end-data.mpsizes.nbl+1:end,[data.linkConst.xu0f.idx data.linkConst.p.idx data.linkConst.t.idx])=1;
data.hessianStruct([data.linkConst.xu0f.idx data.linkConst.p.idx data.linkConst.t.idx],[data.linkConst.xu0f.idx data.linkConst.p.idx data.linkConst.t.idx])=1;
data.hessianStruct=tril(data.hessianStruct);


data.linkConst.all.idx=repmat([data.linkConst.xu0f.idx data.linkConst.p.idx data.linkConst.t.idx],data.mpsizes.nbl,1);


infoNLP_T.zl=problem.mp.time.t_min';
infoNLP_T.zu=problem.mp.time.t_max';
infoNLP_T.z0=guess.mp.time';
mpinfoNLP.zl=[infoNLP_XU.zl;infoNLP_P.zl;infoNLP_T.zl];
mpinfoNLP.zu=[infoNLP_XU.zu;infoNLP_P.zu;infoNLP_T.zu];
mpinfoNLP.z0=[infoNLP_XU.z0;infoNLP_P.z0;infoNLP_T.z0];
mpinfoNLP.cl=[infoNLP_XUg.cl;infoNLP_b.cl;problem.mp.constraints.bll.linear';problem.mp.constraints.bll.nonlinear'];
mpinfoNLP.cu=[infoNLP_XUg.cu;infoNLP_b.cu;problem.mp.constraints.blu.linear';problem.mp.constraints.blu.nonlinear'];

data.multipliers.lambda=[];


auxdata.mpdata=data;
auxdata.mpdata.options=options;

auxdata.mpdata.funcs.hessianstructure  = @hessianstructure_mp;
auxdata.mpdata.funcs.hessian           = @computeHessian_mp;
auxdata.mpdata.funcs.objective         = @costFunction_mp;
auxdata.mpdata.funcs.gradient          = @costGradient_mp;
auxdata.mpdata.funcs.constraints       = @constraintFunction_mp;
auxdata.mpdata.funcs.jacobian          = @constraintJacobian_mp;
auxdata.mpdata.funcs.jacobianstructure = @jacobianstructure_mp;
auxdata.mpdata.linkfunctions=problem.mp.linkfunctions;

auxdata.phasedata=phasedata;
infoNLP.mpinfoNLP=mpinfoNLP;
infoNLP.phaseinfoNLP=phaseinfoNLP;




end

