function [infoNLP,data,options]=transcribeOCP_eachPhase(problem,guess,options)
%TRANSCRIBEOCP - Process information from 'problem', 'guess' and 'options' for NLP solver
%Specifically:
%Error checking of function definitions and bounds
%Define bounds for NLP variable + continuity, path and boundary constraints
%Format matrices for direct transcription method(if required)
%Generate initial guess for optimization
%Generate structure of the jacobian of the constraints(if required)
%Construct optimal finite-difference pertubation sets(if required)
%
% Syntax:  [infoNLP,data,options]=transcribeOCP_eachPhase(problem,guess,options)
%
% Inputs:
%    problem - Optimal control problem definition
%    guess   - Guess used to generate starting point for optimization
%    options - Settings in file settings.m
%
% Outputs:
%    infoNLP - Information required by the NLP solver
%    data - Data passed to the functions evaluated during optimization
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

persistent adigatorGen

checkProblem(problem);

if ~isfield(problem.data,'mode')
    problem.data.mode.currentMode='Original';
end

% Discretization Error Tol
options.discErrorTol=problem.states.xErrorTol_local;
if strcmp(options.transcription,'integral_res_min') || strcmp(options.transcription,'direct_collocation_intres_reg')
    discErrorTol_Full=[problem.states.xErrorTol_integral problem.constraints.gTol_eq];
    global mode_min_res
    mode_min_res=1;
end

% Constraint Violation Tol
options.constraintErrorTol=problem.constraintErrorTol;

% Automatic discretization 
if (strcmp(options.discretization,'AutoDirect')) && ~isfield(options,'AutoDirect')
    if ~strcmp(options.discretization,'manual')
        options.AutoDirect=1;
        options.discretization='hermite';
    else
        error('automatic mesh selection does not support manually chosen solution representation methods')
    end
end

% Adaptively spaced segments
if strcmp(options.meshstrategy,'hp_flexible')
    options.adaptseg=1; 
else
    options.adaptseg=0;
end

if strcmp(options.discretization,'hpLGR') && options.adaptseg
    options.adaptseg=1;
    if isfield(options,'nsegment')
        problem.constraints.bu=[problem.constraints.bl,ones(1,options.nsegment)*options.maxtimeinterval];
        problem.constraints.bl=[problem.constraints.bl,ones(1,options.nsegment)*options.mintimeinterval];
    elseif isfield(options,'npsegment')
        problem.constraints.bu=[problem.constraints.bl,ones(1,length(options.npsegment))*options.maxtimeinterval];
        problem.constraints.bl=[problem.constraints.bl,ones(1,length(options.npsegment))*options.mintimeinterval];
    end
end

% Define (and assign) some parameters
%---------------------------------------
n=length(problem.states.x0l);              % Number of states
m=length(problem.inputs.ul);               % Number of inputs
np=length(problem.parameters.pl);          % Number of free parameters
ng_neq=length(problem.constraints.gl);         % Number of path constraints
ng_eq=problem.constraints.ng_eq;         % Number of path constraints
nb=length(problem.constraints.bl);         % Number of boundary constraints
N=problem.inputs.N;                        % Number of control actions



if ~isfield(problem.states,'resNormCusWeight') && ~isfield(problem.constraints,'resNormCusWeight_eq')
    problem.data.resNormCusWeight=ones(1,n+ng_eq);
elseif isfield(problem.states,'resNormCusWeight') && isfield(problem.constraints,'resNormCusWeight_eq')
    problem.data.resNormCusWeight=[problem.states.resNormCusWeight.*ones(1,n) problem.constraints.resNormCusWeight_eq.*ones(1,ng_eq)];
elseif ~isfield(problem.states,'resNormCusWeight') && isfield(problem.constraints,'resNormCusWeight_eq')
    problem.data.resNormCusWeight=[ones(1,n) problem.constraints,resNormCusWeight_eq.*ones(1,ng_eq)];
elseif isfield(problem.states,'resNormCusWeight') && ~isfield(problem.constraints,'resNormCusWeight_eq')
    problem.data.resNormCusWeight=[problem.states.resNormCusWeight.*ones(1,n) ones(1,ng_eq)];
else
    error('Customized weighting for residual norm not properly configured! Please assign the relative weighting with problem.states.resNormCusWeight and problem.constraints,resNormCusWeight_eq accordingly')
end

data.data=problem.data;
if ~isfield(problem, 'analyticDeriv')
    problem.analyticDeriv.gradCost=[];
    problem.analyticDeriv.hessianLagrangian=[];
    problem.analyticDeriv.jacConst=[];
end
data.analyticDeriv = problem.analyticDeriv;
data.data.transcription=options.transcription;
data.data.discretization=options.discretization;

% Scale problem variables
if options.scaling 
    [ problem,guess ] = scale_problem( problem,guess );
end

% Variable rate constraints
if isfield(problem.states,'xrl')
    data.data.xrl=problem.states.xrl;
    data.data.xru=problem.states.xru;
    data.data.xrConstraintTol=problem.states.xrConstraintTol;
    data.data.url=problem.inputs.url;
    data.data.uru=problem.inputs.uru;
    data.data.urConstraintTol=problem.inputs.urConstraintTol;
end

% Two options for the LGR points on interval t=[-1,tau_n], tau_n<1. The
% corresponding parameters are extracted. 
% variables: 
%   npd:        a vector contains the polynomial degree of each segment
%   npdu:       unique entries of npd
%   npduidx:    a vector same length as npd, indicating the corresponding
%               indeix in npdu
%   nps:        number of segments (scalar)
%   M:          total number of mesh nodes (scalar)
if (strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))
    if isfield(options,'pdegree') 
        npd=options.pdegree*ones(1,options.nsegment);
        [npdu,~,npduidx]=unique(npd);
        nps=options.nsegment;
        M=sum(npd); 
    else
        npd=options.npsegment;
        [npdu,~,npduidx]=unique(npd);                  
        nps=length(options.npsegment);                
        M=sum(npd);                            
    end
    config_Check(options,nps,npd);
else
    if length(options.tau)~=1
        M=length(options.tau)+1;
    else    
        M=options.nodes;  
    end
    nps=M-1;
    config_Check(options,M);
end


% Multiple-shooting: Number of integration nodes in the interval t=[t0,tf] is N+1
if strcmp(options.transcription,'multiple_shooting')
   if N==0
       N=M-1;
   else
       M=N+1;
   end % Default N=M-1.
else
    if isfield(options,'disContInputs')&& options.disContInputs
        if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
            if N==0
                N=M-1;
                ns=1;
            end
        elseif strcmp(options.discretization,'hermite')
            if N==0
                N=(M-1)*3;
                M=2*M-1;
                ns=2;
            end
        elseif strcmp(options.discretization,'trapezoidal')
            if N==0
                N=M*2-2;
                ns=1;
            end
        elseif strcmp(options.discretization,'euler')
            if N==0
                N=M;
                ns=1;
            end
        end
        
    else
        if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
            if N==0
                N=M-1;
                ns=1;
            end
        elseif strcmp(options.discretization,'hermite')
            if N==0
                N=M*2-1;
                M=2*M-1;
                ns=2;
            end
        elseif strcmp(options.discretization,'discrete')
            if N==0
                N=M;
                ns=1;
            end
        else
            if N==0
                N=M;
                ns=1;
            end
        end
        
    end

end


% The final time for discrete time systems is imposed equal to 1 in order 
% to use a unform formulation of the optimization problem
if (strcmp(options.discretization,'discrete'))
   problem.time.tf_min=1;
   problem.time.tf_max=1;
   guess.tf=1;
end 
 

% Get bounds for final time and check if time is free or fixed
% nt=0 when the final time is not a variable otherwise nt=1.
% For LGR, when t0 and tf both variables: nt=2
%          when adaptive method, nt equal to the actual number of time
%          varibles (>2)

if (strcmp(options.NLPsolver,'OSQP')) 
   if problem.time.tf_min==problem.time.tf_max && problem.time.t0_min==problem.time.t0_max
       options.runWithoutTimeVar=1;
   else
       error('To solve with OSQP, the problem must have fixed intial time t0 and terminal time tf.')
   end
end 

if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
    if options.adaptseg==1
        nt=nps+1;
    else
        nt=2;
    end
else
    if isfield(options,'runWithoutTimeVar') && options.runWithoutTimeVar
        nt=0;
    else
        nt=2; 
    end
end

if nt
    t0l=problem.time.t0_min; t0u=problem.time.t0_max;
    tfl=problem.time.tf_min; tfu=problem.time.tf_max;
else
    t0l=[];t0u=[];tfl=[];tfu=[];
end

% Set the initial time and other parameters to adjust the temporal
% scale

if ~nt
    data.t0=problem.time.t0_min;
    data.tf=problem.time.tf_min;
end

if (strcmp(options.discretization,'discrete'))
     if problem.time.t0_min~=problem.time.t0_max
        error('Error: for discrete problems, currently only support fixed t0. Please set t0_min to equal to t0_max to precede.')
     else
         data.k0=problem.time.t0_min;
         problem.time.t0=0;
         data.Nm=N;
     end
else
     data.Nm=1;
end



% Number of algebraic variable rate constraints
if (strcmp(options.discretization,'hermite'))
    if isfield(problem.states,'xrl')
        
        rcxl=find(~isinf(problem.states.xrl) .* (problem.states.xrl~=problem.states.xru));
        rcxu=find(~isinf(problem.states.xru) .* (problem.states.xrl~=problem.states.xru));
        rcul=find(~isinf(problem.inputs.url) .* (problem.inputs.url~=problem.inputs.uru));
        rcuu=find(~isinf(problem.inputs.uru) .* (problem.inputs.url~=problem.inputs.uru));
        rcxe=find(problem.states.xrl==problem.states.xru);
        rcue=find(problem.inputs.url==problem.inputs.uru);

        data.RCmap.rcxl=rcxl;
        data.RCmap.rcxu=rcxu;
        data.RCmap.rcul=rcul;
        data.RCmap.rcuu=rcuu;
        data.RCmap.rcxe=rcxe;
        data.RCmap.rcue=rcue;
        
        vnrcl=(sum(~isinf(problem.states.xrl))-sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru))))*3+(sum(~isinf(problem.inputs.url))-sum(problem.inputs.url(~isinf(problem.inputs.url))==problem.inputs.uru(~isinf(problem.inputs.uru))))*2;
        vnrcu=(sum(~isinf(problem.states.xru))-sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru))))*3+(sum(~isinf(problem.inputs.uru))-sum(problem.inputs.url(~isinf(problem.inputs.url))==problem.inputs.uru(~isinf(problem.inputs.uru))))*2;
        vnrce=sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru)))*2+sum(problem.inputs.url(~isinf(problem.inputs.url))==problem.inputs.uru(~isinf(problem.inputs.uru)))*2;
        nrcl=(sum(~isinf(problem.states.xrl))-sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru))))*3*(nps)+(sum(~isinf(problem.inputs.url))-sum(problem.inputs.url(~isinf(problem.inputs.url))==problem.inputs.uru(~isinf(problem.inputs.uru))))*2*(nps);
        nrcu=(sum(~isinf(problem.states.xru))-sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru))))*3*(nps)+(sum(~isinf(problem.inputs.uru))-sum(problem.inputs.url(~isinf(problem.inputs.url))==problem.inputs.uru(~isinf(problem.inputs.uru))))*2*(nps);
        nrce=sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru)))*2*(nps)+sum(problem.inputs.url(~isinf(problem.inputs.url))==problem.inputs.uru(~isinf(problem.inputs.uru)))*2*(nps);
    else
        vnrcl=0;vnrcu=0;vnrce=0;nrcl=0;nrcu=0;nrce=0;
    end

elseif strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
    if isfield(problem.states,'xrl')
        rcxl=find(~isinf(problem.states.xrl) .* (problem.states.xrl~=problem.states.xru));
        rcxu=find(~isinf(problem.states.xru) .* (problem.states.xrl~=problem.states.xru));
        rcul=find(~isinf(problem.inputs.url) .* (problem.inputs.url~=problem.inputs.uru));
        rcuu=find(~isinf(problem.inputs.uru) .* (problem.inputs.url~=problem.inputs.uru));
        rcxe=find(problem.states.xrl==problem.states.xru);
        rcue=find(problem.inputs.url==problem.inputs.uru);
        
        vnrcl=length(rcxl)+length(rcul);
        vnrcu=length(rcxu)+length(rcuu);
        vnrce=length(rcxe)+length(rcue);
        nrcl=length(rcxl)*M+length(rcul)*M;
        nrcu=length(rcxu)*M+length(rcuu)*M;
        nrce=length(rcxe)*M+length(rcue)*M;
        
        data.RCmap.rcxl=rcxl;
        data.RCmap.rcxu=rcxu;
        data.RCmap.rcul=rcul;
        data.RCmap.rcuu=rcuu;
        data.RCmap.rcxe=rcxe;
        data.RCmap.rcue=rcue;
    else
        vnrcl=0;vnrcu=0;vnrce=0;nrcl=0;nrcu=0;nrce=0;
    end
else 
    if isfield(problem.states,'xrl')
        rcxl=find(~isinf(problem.states.xrl) .* (problem.states.xrl~=problem.states.xru));
        rcxu=find(~isinf(problem.states.xru) .* (problem.states.xrl~=problem.states.xru));
        rcul=find(~isinf(problem.inputs.url) .* (problem.inputs.url~=problem.inputs.uru));
        rcuu=find(~isinf(problem.inputs.uru) .* (problem.inputs.url~=problem.inputs.uru));
        rcxe=find(problem.states.xrl==problem.states.xru);
        rcue=find(problem.inputs.url==problem.inputs.uru);
        
        data.RCmap.rcxl=rcxl;
        data.RCmap.rcxu=rcxu;
        data.RCmap.rcul=rcul;
        data.RCmap.rcuu=rcuu;
        data.RCmap.rcxe=rcxe;
        data.RCmap.rcue=rcue;
        
        vnrcl=length(rcxl)+length(rcul);
        vnrcu=length(rcxu)+length(rcuu);
        vnrce=length(rcxe)+length(rcue);
        nrcl=length(rcxl)*(M-1)+length(rcul)*(N-1);
        nrcu=length(rcxu)*(M-1)+length(rcuu)*(N-1);
        nrce=length(rcxe)*(M-1)+length(rcue)*(N-1);
    else
        vnrcl=0;vnrcu=0;vnrce=0;nrcl=0;nrcu=0;nrce=0;
    end
end
nrc=nrcl+nrcu+nrce; 


if nrc
    if isfield(options,'disContInputs') && options.disContInputs
       error('Rate constraints not supported when having discontinuous inputs') 
    end
    
    if strcmp(options.discretization,'hermite')

        data.RCmap.AxHS1=spdiags([-3*ones(M,1) 4*ones(M,1) -ones(M,1)],0:2,M,M);
        data.RCmap.AxHS1([2:2:end,end],:)=[];
        data.RCmap.AxHS2=spdiags([ones(M,1) -4*ones(M,1) 3*ones(M,1)],0:2,M,M);
        data.RCmap.AxHS2([2:2:end,end],:)=[];
        data.RCmap.AxHS3=spdiags([-ones(M,1) ones(M,1)],[0 2],M,M);
        data.RCmap.AxHS3([2:2:end,end],:)=[];
        data.RCmap.AxHS4=spdiags([-ones(M,1) ones(M,1)],[0 1],M,M);
        data.RCmap.AxHS4([2:2:end,end],:)=[];

        data.RCmap.AuHS1=spdiags([-3*ones(N,1) 4*ones(N,1) -ones(N,1)],0:2,N,N);
        data.RCmap.AuHS1([2:2:end,end],:)=[];
        data.RCmap.AuHS2=spdiags([ones(N,1) -4*ones(N,1) 3*ones(N,1)],0:2,N,N);
        data.RCmap.AuHS2([2:2:end,end],:)=[];
        data.RCmap.AuHS3=spdiags([-ones(N,1) ones(N,1)],[0 2],N,N);
        data.RCmap.AuHS3([2:2:end,end],:)=[];
        data.RCmap.AuHS4=spdiags([-ones(N,1) ones(N,1)],[0 1],N,M);
        data.RCmap.AuHS4([2:2:end,end],:)=[];
    elseif (strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))
    elseif strcmp(options.discretization,'trapezoidal')
        data.RCmap.Ax=spdiags([-ones(M,1) ones(M,1)],[0 1],M,M);
        data.RCmap.Ax(end,:)=[];
        data.RCmap.Au=spdiags([-ones(N,1) ones(N,1)],[0 1],N,M);
        data.RCmap.Au(end,:)=[];
    else
        error('Rate constraints not supported with the chosen discretization method') 
    end
end

% Generation of mesh (time dimension)
if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
    [ data, tau_inc, tau_seg, tau, LGR ] = genTimeMeshLGR( problem,options, data, nps, npdu, npd, npduidx, M );
else
    if strcmp(options.discretization,'discrete')
        [ data, tau ] = genTimeMesh( problem, options, data, ns, M );
    else
        [ data, tau ] = genTimeMesh( problem, options, data, ns, M );
    end
end



if ng_neq
    if isfield(problem.constraints,'g_neq_ActiveTime') && problem.time.t0_min==problem.time.t0_max && problem.time.tf_min==problem.time.t0_max %%&& options.ECH.enabled 
        if (strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')) && options.reorderLGR && ~all(cellfun('isempty',problem.constraints.g_neq_ActiveTime))
            error('Alternative ordering of variables with LGR method currently does not support external constraint handling in time dimension')
        else
            for i=1:ng_neq
                if ~isempty(problem.constraints.g_neq_ActiveTime{i})
                    gActiveTime=problem.constraints.g_neq_ActiveTime{i}/(guess.tf-guess.t0)*ns;
                    for j=1:size(problem.constraints.g_neq_ActiveTime{i},1)
                        if (strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))
                            gActiveTimeLower=data.tau_segment(data.tau_segment<=(gActiveTime(j,1)*2-1));
                            gActiveTimeUpper=data.tau_segment(data.tau_segment>=(gActiveTime(j,2)*2-1));
                            gActiveIdxi(:,j)=(data.tau_inc>=gActiveTimeLower(end) & data.tau_inc<=gActiveTimeUpper(1));
                        else
                            gActiveIdxi(:,j)=([0;data.tau_inc]>=gActiveTime(j,1) & [0;data.tau_inc]<=gActiveTime(j,2));
                        end
                    end
                    gActiveIdx(:,i)=logical(sum(gActiveIdxi,2));
                else
                    gActiveIdx(:,i)=false(M,1);
                end
            end
            ngActive=sum(sum(gActiveIdx));
        end
    else
        gActiveIdx=true(M,ng_neq);
        ngActive=M*ng_neq;
    end
else
    gActiveIdx=[];
    ngActive=0;
end

if strcmp(options.transcription,'integral_res_min') 
    ng=ng_neq;
    problem.constraints.gTol=[problem.constraints.gTol_neq];
else
    ng=ng_eq+ng_neq;
    gActiveIdx=[true(M,ng_eq) gActiveIdx];
    ngActive=ngActive+M*ng_eq;
    problem.constraints.gl=[zeros(1,ng_eq) problem.constraints.gl];
    problem.constraints.gu=[zeros(1,ng_eq) problem.constraints.gu];
    problem.constraints.gTol=[problem.constraints.gTol_eq problem.constraints.gTol_neq];
end
            
if isfield(problem.data,'gFilter')
    data.data.gFilter=problem.data.gFilter;
end

% Number of path constraints
% Save the corrsponding variables sizes
if (strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))
    data.sizes={nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive,ng_eq,ng_neq};
else    
    data.sizes={nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive,nps,ng_eq,ng_neq};
end

% Get state and input bounds
xl=problem.states.xl(:); xu=problem.states.xu(:);
ul=problem.inputs.ul(:); uu=problem.inputs.uu(:);

% Choose the smallest constraint set for the terminal bounds
xfl=problem.states.xfl(:);xfu=problem.states.xfu(:);
xf_l=xl;xf_l(xl<xfl)=xfl(xl<xfl);
xf_u=xu;xf_u(xu>xfu)=xfu(xu>xfu);

% Choose the smallest constraint set for the initial state bounds
x0l=problem.states.x0l(:);x0u=problem.states.x0u(:);
x0_l=xl;
x0_l(xl<x0l)=x0l(xl<x0l);
x0_u=xu;
x0_u(xu>x0u)=x0u(xu>x0u);
u0_l=problem.inputs.u0l(:);u0_u=problem.inputs.u0u(:);

% configure the initial condition
if ~isfield(problem.states,'x0') || isempty(problem.states.x0)
  if isfield(guess,'states')
      cx0=0;
      data.x0=guess.states(1,:);
      data.x0t=guess.states(1,:)';
  else
      cx0=0;
      data.x0=x0l;
      data.x0t=(x0_l+(x0_u-x0_l).*rand(n,1));
      data.x0t(isinf(data.x0t))=0;
      data.x0t(isnan(data.x0t))=0;
  end
else  
  cx0=1;
  data.x0t=problem.states.x0.';
  data.x0=problem.states.x0;
end 
data.cx0=cx0;


% Get parameter bounds
if np
    pl=problem.parameters.pl(:); 
    pu=problem.parameters.pu(:);
else 
    pl=[];
    pu=[]; 
end

% Get bounds for the constraint functions
if ng 
    gl=problem.constraints.gl(:);
    gu=problem.constraints.gu(:);
    g_tol=problem.constraints.gTol(:);
    if (strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')) && length(gl)>1
        glAll=repelem(gl,M);
        guAll=repelem(gu,M);
        g_tolAll=repelem(g_tol,M);
        if ng_eq
            g_tolEq=problem.constraints.gTol_eq(:);
            g_tolEq=repelem(g_tolEq,M);
        else
            g_tolEq=[];
        end
        gAllidx=logical(gActiveIdx(:));
    elseif strcmp(options.transcription,'multiple_shooting')
        glAll=kron(ones(M-1,1),gl(:));
        guAll=kron(ones(M-1,1),gu(:));
        g_tolAll=kron(ones(M-1,1),g_tol(:));
%         gAllidx=logical(reshape(gActiveIdx',size(gActiveIdx,1)*size(gActiveIdx,2),1));
    else
        
        glAll=kron(ones(M,1),gl(:));
        guAll=kron(ones(M,1),gu(:));
        g_tolAll=kron(ones(M,1),g_tol(:));
        if ng_eq
            g_tolEq=problem.constraints.gTol_eq(:);
            g_tolEq=kron(ones(M,1),g_tolEq(:));
        else
            g_tolEq=[];
        end
        gAllidx=logical(reshape(gActiveIdx',size(gActiveIdx,1)*size(gActiveIdx,2),1));
    end
    data.gActiveIdx=gActiveIdx;
else 
    gl=[];
    gu=[];
    glAll=[];
    guAll=[];
    gAllidx=[];
    g_tolAll=[];
    g_tolEq=[];
    data.gActiveIdx=gActiveIdx;

end
data.glAll=glAll;
data.guAll=guAll;
data.g_tolAll=g_tolAll;
data.g_tolEq=g_tolEq;
data.gAllidx=find(gAllidx);



% Get bounds for the constraint functions
if nrc
    rcl=[zeros(1,nrcl),-inf*ones(1,nrcu),zeros(1,nrce)];
    rcu=[inf*ones(1,nrcl),zeros(1,nrcu),zeros(1,nrce)];
else 
    rcl=[];
    rcu=[];
end

if nb
    bl=problem.constraints.bl(:);
    bu=problem.constraints.bu(:);
    bTol=problem.constraints.bTol(:);
else 
    bl=[];
    bu=[];
    bTol=[];
end

        
% Store some matrices in data structure
data.options=options;
data.functions=problem.functions;
data.functions_unscaled=problem.functions_unscaled;
data.data.N_tNode=M;

if options.scaling
    data.data.Xscale=problem.states.scales;
    data.data.Uscale=problem.inputs.scales;
    data.data.Xscale_back=problem.states.scales_back;
    data.data.Uscale_back=problem.inputs.scales_back;
    data.data.Tscale=problem.time.scales;
    data.data.Tshift=problem.time.shifts;
    data.data.Xshift=problem.states.shifts;
    data.data.Ushift=problem.inputs.shifts;
    
    
    if np
        data.data.Pscale=problem.parameters.scales;
        data.data.Pscale_back=problem.parameters.scales_back;
        data.data.Pshift=problem.parameters.shifts;
        if isfield(options,'runWithoutTimeVar') && options.runWithoutTimeVar && nt==0
            data.data.Allscale_fgL=[data.data.Xscale data.data.Uscale data.data.Pscale ];
            data.data.Allscale_bE=[data.data.Xscale data.data.Xscale data.data.Uscale data.data.Uscale data.data.Pscale ];
        else
            data.data.Allscale_fgL=[data.data.Xscale data.data.Uscale data.data.Pscale data.data.Tscale data.data.Tscale];
            data.data.Allscale_bE=[data.data.Xscale data.data.Xscale data.data.Uscale data.data.Uscale data.data.Pscale data.data.Tscale data.data.Tscale];
        end
    else
        if isfield(options,'runWithoutTimeVar') && options.runWithoutTimeVar && nt==0
            data.data.Allscale_fgL=[data.data.Xscale data.data.Uscale ];
            data.data.Allscale_bE=[data.data.Xscale data.data.Xscale data.data.Uscale data.data.Uscale ];
        else
            data.data.Allscale_fgL=[data.data.Xscale data.data.Uscale data.data.Tscale data.data.Tscale];
            data.data.Allscale_bE=[data.data.Xscale data.data.Xscale data.data.Uscale data.data.Uscale data.data.Tscale data.data.Tscale];
        end

    end
    data.data.Allscale_fgL_Mat=data.data.Allscale_fgL'*data.data.Allscale_fgL;
    data.data.Allscale_bE_Mat=data.data.Allscale_bE'*data.data.Allscale_bE;
    
    if strcmp(options.derivatives,'adigator') || strcmp(options.transcription,'integral_res_min') || strcmp(options.discretization,'resMinInterpolationForSolution') 
        if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
            XunscaleMat=1./repmat( data.data.Xscale, M+1, 1 );
            data.scaling.XunscaleMat=spdiags(XunscaleMat(:),0,(M+1)*n,(M+1)*n);
            XscaleMat=repmat( data.data.Xscale, (M+1), 1 );
            data.scaling.XscaleMat=spdiags(XscaleMat(:),0,(M+1)*n,(M+1)*n);
            data.scaling.XshiftMat=repmat( data.data.Xshift, (M+1), 1 );

            UunscaleMat=1./repmat( data.data.Uscale, M, 1 );
            data.scaling.UunscaleMat=spdiags(UunscaleMat(:),0,M*m,M*m);
            UscaleMat=repmat( data.data.Uscale, M, 1 );
            data.scaling.UscaleMat=spdiags(UscaleMat(:),0,M*m,M*m);
            data.scaling.UshiftMat=repmat( data.data.Ushift, M, 1 );
            if isfield(data.data,'Pscale')
                PscaleMat=repmat( data.data.Pscale, M, 1 );
                PunscaleMat=1./repmat( data.data.Pscale, M, 1 );
                data.scaling.PshiftMat=repmat( data.data.Pshift, M, 1 );
                data.scaling.PunscaleMat=spdiags(PunscaleMat(:),0,M*np,M*np);
                data.scaling.PscaleMat=spdiags(PscaleMat(:),0,M*np,M*np);
            end
        else
            XunscaleMat=1./repmat( data.data.Xscale, M, 1 );
            data.scaling.XunscaleMat=spdiags(XunscaleMat(:),0,M*n,M*n);
            XscaleMat=repmat( data.data.Xscale, M, 1 );
            data.scaling.XscaleMat=spdiags(XscaleMat(:),0,M*n,M*n);
            data.scaling.XshiftMat=repmat( data.data.Xshift, M, 1 );

            UunscaleMat=1./repmat( data.data.Uscale, M, 1 );
            data.scaling.UunscaleMat=spdiags(UunscaleMat(:),0,M*m,M*m);
            UscaleMat=repmat( data.data.Uscale, M, 1 );
            data.scaling.UscaleMat=spdiags(UscaleMat(:),0,M*m,M*m);
            data.scaling.UshiftMat=repmat( data.data.Ushift, M, 1 );

            if isfield(data.data,'Pscale')
                PscaleMat=repmat( data.data.Pscale, M, 1 );
                PunscaleMat=1./repmat( data.data.Pscale, M, 1 );
                data.scaling.PshiftMat=repmat( data.data.Pshift, M, 1 );
                data.scaling.PunscaleMat=spdiags(PunscaleMat(:),0,M*np,M*np);
                data.scaling.PscaleMat=spdiags(PscaleMat(:),0,M*np,M*np);
            end
        end
    end
    

end
if isempty(data.options.perturbation.J)
  data.options.perturbation.J=(eps/2)^(1/3);
end
if isempty(data.options.perturbation.H)    
  data.options.perturbation.H=(8*eps)^(1/3);
end    

% Reference trajectory for the cost function; it is assigned as the set-point
% but the user can be redefine it if necessary
if isfield(problem,'setpoints')
    data.references.xr=repmat(problem.setpoints.states,M,1);   
    data.references.ur=repmat(problem.setpoints.inputs,M,1);   
else
    data.references.xr=[];
    data.references.ur=[];
end

% Define bounds for the NLP variable 
%---------------------------------------
% A) direct multiple shooting              : z = [tf p x(0) u(0) x(1) u(1) ... x(M)]'
%
% B) h method direct transcription (N = M-1) :
%    z = [tf p x(0) u(0) x(1) u(1) ... x(M-1) u(M-1) x(M) ]' 
%    In the trapezoidal method the required u(M) is imposed equal to u(M-1)
%
% C) h method direct transcription (N < M-1) : z = [tf p x(0) x(1) ...x((M-1)/N-1) u(0) x((M-1)/N)...x(M-1) u(N-1) x(M)]'
% D) p(LGR) method direct transcription: z = [x(1) x(2) ... x(M+1) u(1) u(2) ... u(M) p t0 ... tf]' 
% E) p(LGR) method direct transcription (alternative): z = [x(1) u(1) x(2) u(2) ... x(M) u(M) x(M+1) p t0 ... tf]' 
%    Enable this by setting options.reorderLGR to 1
switch options.transcription
        
    case{'multiple_shooting'}
        if N<1; error('Number of control actions incorrect'); end
        nx=n*M;                         % Number of unknown states
        nu=N*m;                         % Number of unknown controls
        nz=nt+np+nx+nu;                 % Length of the optimization variable
        infoNLP.zl=[tfl; pl;  kron(ones(N,1),[xl(:);ul(:)]);xf_l(:)];
        infoNLP.zu=[tfu; pu;  kron(ones(N,1),[xu(:);uu(:)]);xf_u(:)];
        infoNLP.zl(nt+np+1:nt+np+n)=x0_l;    % Prune bounds for initial conditions
        infoNLP.zu(nt+np+1:nt+np+n)=x0_u;
        data.Nm=1;
        data.options.ipopt.hessian_approximation='limited-memory';
        
        data.zidx.org.xu=nt+np+1:nz;
        data.zidx.org.t=1:nt;
        data.zidx.org.p=nt+1:nt+np;
    
    case{'direct_collocation','direct_collocation_intres_reg','integral_res_min'}
        
        switch options.discretization
            case{'discrete','euler','trapezoidal','hermite'}
                if mod(M,N)||(N>M)              % Check for errors
                    error('# integration steps is not divisible or less than N');
                end
        %         nk=M-1;                      % Number of integration steps
                nx=M*n;                      % Number of unknown states  
                nu=N*m;                      % Number of unknown controls
                nz=nt+np+nx+nu;              % Length of the optimization variable

                xpl=repmat(xl(:),M/N,1);xpu=repmat(xu(:),M/N,1);
                
                infoNLP.zl=[t0l;tfl; pl;repmat([xpl(:);ul(:)],N-1,1);xf_l(:);ul(:)];
                infoNLP.zu=[t0u;tfu; pu;repmat([xpu(:);uu(:)],N-1,1);xf_u(:);uu(:)];
                infoNLP.zl(nt+np+1:nt+np+n)=x0_l;    % Prune bounds for initial conditions
                infoNLP.zu(nt+np+1:nt+np+n)=x0_u;
                infoNLP.zl(nt+np+n+1:nt+np+n+m)=u0_l;    % Prune bounds for initial conditions
                infoNLP.zu(nt+np+n+1:nt+np+n+m)=u0_u;
                
                data.zidx.org.xu=nt+np+1:nz;
                data.zidx.org.t=1:nt;
                data.zidx.org.p=nt+1:nt+np;

            case{'globalLGR','hpLGR'}
                if mod(M-1,N)||(N>M)              % Check for errors
                    error('# integration steps + 1 not divisible or less than N');
                end
        %         nk=M;                        % Number of integration steps
                nx=(M+1)*n;                  % Number of unknown states  
                nu=M*m;                      % Number of unknown controls
                nz=nt+np+nx+nu;              % Length of the optimization variable

                xpl=repmat(xl(:),(M-1)/N,1);xpu=repmat(xu(:),(M-1)/N,1);

                xl_mat=[x0_l';repmat(xpl(:)',M-1,1);xf_l']; %lower bounds of X
                ul_mat=[u0_l';repmat(ul(:)',M-1,1);]; %lower bounds of U
                xu_mat=[x0_u';repmat(xpu(:)',M-1,1);xf_u']; %upper bounds of X
                uu_mat=[u0_u';repmat(uu(:)',M-1,1);]; %upper bounds of U
                if options.adaptseg==1
                    tl_mat=[t0l cumsum(ones(1,nt-1)*options.mintimeinterval)];
                    tu_mat=[t0u cumsum(ones(1,nt-1)*options.maxtimeinterval)];
                else
                    tl_mat=[t0l tfl];
                    tu_mat=[t0u tfu];
                end       
                infoNLP.zl=[xl_mat(:);ul_mat(:);pl;tl_mat']; %lower bounds of z
                infoNLP.zu=[xu_mat(:);uu_mat(:);pu;tu_mat']; %upper bounds of z
                
                data.zidx.org.xu=1:nx+nu;
                data.zidx.org.t=nz-nt+1:nz;
                data.zidx.org.p=nz-nt-np+1:nz-nt;
            otherwise;error('Unknown discretization method. Check spelling');
        end
    otherwise;error('Unknown transcription method. Check spelling');
end
data.nz=nz;

if strcmp(options.transcription,'integral_res_min') && (strcmp(options.discretization,'discrete') || strcmp(options.discretization,'euler') || strcmp(options.discretization,'trapezoidal'))
    error('Selected discretization method not supported by integral residual minimization. Please use direct collocation as transcription method, or select a different discretization method');
end



% Define bounds for constraint functions
%---------------------------------------
% Constraints: c = [c0... cF g(0)....g(ng-1) bo, ,bo(nb-1)]' 
%              ck -> defects  for k=0, ..., F=M-1 
%              gk -> general path constraints g(x(k),u(k),p,k) for
%              k=0,...,ng-1
%              bo(k) -> nb boundary conditions b(x(0),x(f),u(0),u(f),p,t)
if isfield(options,'directColleqtol')
    directColleqtol=options.directColleqtol;
else
    directColleqtol=eps;
end
% directColleqtol=1e-9;
if strcmp(options.transcription,'multiple_shooting')
    infoNLP.cl=[kron(ones(M,1),zeros(n,1));glAll;bl(:)];
    infoNLP.cu=[kron(ones(M,1),zeros(n,1));guAll;bu(:)];   
elseif strcmp(options.transcription,'direct_collocation') || strcmp(options.transcription,'direct_collocation_intres_reg')
    infoNLP.cl=[kron(ones(M,1),-directColleqtol*ones(n,1));glAll(gAllidx);rcl(:);bl(:)];
    infoNLP.cu=[kron(ones(M,1),directColleqtol*ones(n,1));guAll(gAllidx);rcu(:);bu(:)];
    if strcmp(options.discretization,'discrete') || strcmp(options.discretization,'euler')
        infoNLP.cl(end-length([rcl(:);bl(:)])-ng+1:end-length([rcl(:);bl(:)]),1)=0;
        infoNLP.cu(end-length([rcl(:);bl(:)])-ng+1:end-length([rcl(:);bl(:)]),1)=0;
    end
    if strcmp(options.transcription,'direct_collocation_intres_reg')
        discErrorTol_Full=discErrorTol_Full.^2;
        discErrorTol_Full=discErrorTol_Full';
        data.data.discErrorTol_Full=discErrorTol_Full;
        data.data.discErrorConstScaling=ones(1,n+ng_eq);
        data.data.discErrorTol_FullScaling=discErrorTol_Full;
    end
%     infoNLP.cl=[kron(ones(M,1),zeros(n,1));glAll(gAllidx);rcl(:);bl(:)];
%     infoNLP.cu=[kron(ones(M,1),zeros(n,1));guAll(gAllidx);rcu(:);bu(:)];
else
%     if options.scaling 
%         discErrorTol_Full= scale_variables( discErrorTol_Full, data.data.Xscale, 0 ).^2;
%         if isfield(guess,'residual')
%             guess.residual= transpose(scale_variables( guess.residual', data.data.Xscale, 0 )).^2;
%         end
%     else
        discErrorTol_Full=discErrorTol_Full.^2;
%     end
    if any(discErrorTol_Full<eps^2)
        error('Integral of the residual errors squared allowed must be strictly larger than machine precision');
    else
        if strcmp(options.min_res_mode,'directCostMin')
            discErrorTol_Full=discErrorTol_Full';
            data.data.discErrorTol_Full=discErrorTol_Full;
            data.data.discErrorConstScaling=1./sqrt(discErrorTol_Full)';
            data.data.discErrorTol_FullScaling=data.data.discErrorTol_Full.*data.data.discErrorConstScaling';
%             ResConstScaleMat=repmat(data.dataNLP.data.discErrorConstScaling, data.nps, 1 );
%             data.ResConstScaleMat=diag(ResConstScaleMat(:));
            
            infoNLP.cl=[glAll(gAllidx);rcl(:);bl(:);zeros(n+ng_eq,1)];
            infoNLP.cu=[guAll(gAllidx);rcu(:);bu(:);ones(n+ng_eq,1).*sqrt(discErrorTol_Full)];
            infoNLP.c_Tol=[g_tolAll(gAllidx);zeros(size(rcu(:)));bTol(:);zeros(n+ng_eq,1)];

        else
            if isfield(guess,'cost')
                cost_ub=guess.cost.ub;
                cost_lb=guess.cost.lb;
            else
                cost_lb=-inf;
                cost_ub=inf;
            end
            discErrorTol_Full=discErrorTol_Full';
            data.data.discErrorTol_Full=discErrorTol_Full;
            if isfield(guess,'residual')
                discErrorConst=max(discErrorTol_Full/2,guess.residual);
                data.data.discErrorConstScaling=1./sqrt(discErrorConst(:))';
                data.data.discErrorTol_FullScaling=discErrorTol_Full.*data.data.discErrorConstScaling;
                infoNLP.cl=[glAll(gAllidx);rcl(:);bl(:);cost_lb;zeros(n+ng_eq,1)];
                infoNLP.cu=[guAll(gAllidx);rcu(:);bu(:);cost_ub;ones(n+ng_eq,1).*sqrt(discErrorConst)];
            else
                discErrorConst=inf*ones(n+ng_eq,1);
                data.data.discErrorConstScaling=ones(1,n+ng_eq);
                data.data.discErrorTol_FullScaling=discErrorTol_Full;
                infoNLP.cl=[glAll(gAllidx);rcl(:);bl(:);cost_lb;zeros(n+ng_eq,1)];
                infoNLP.cu=[guAll(gAllidx);rcu(:);bu(:);cost_ub;discErrorConst(:)];
            end
            infoNLP.c_Tol=[g_tolAll(gAllidx);zeros(size(rcu(:)));bTol(:);0;zeros(n+ng_eq,1)];
        end
        
    end
end

data.nConst=length(infoNLP.cl);

% Extract sparsity structures
%---------------------------------------
% if ~strcmp(options.derivatives,'analytic')
% 
%     ~(isfield(options,'sysStructTest') && ~options.sysStructTest)
sparsity_num=getStructure(data.sizes,options.discretization);
if ~strcmp(options.derivatives,'adigator') && (strcmp(options.resultRep,'default') || strcmp(options.resultRep,'manual')) && (strcmp(options.transcription,'direct_collocation') || strcmp(options.transcription,'direct_collocation_intres_reg')) && (isfield(options,'sysStructTest') && options.sysStructTest) && ~(strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR'))
    try
        sparsity=getStructure_NaNTest( problem, options, data, sparsity_num );
    catch
        sparsity=sparsity_num;
    end
else
    sparsity=sparsity_num;
end
if strcmp(options.derivatives,'analytic')
    [sparsity, data]=getStructureA(sparsity,data);
end

% else
    % The structure of the derivatives is determined only considering a
    % a fixed structure
    
%     sparsity_num=getStructure(data.sizes,options.discretization);
% end
% data.sparsity=sparsity; 




% Format direct transcription matrices
%---------------------------------------
% Format matrices for continuity constraints: c(z)=[A.Vx.z+B.F(z)]
% and generate the quadrature vector w(tau)

if ~strcmp(options.transcription,'multiple_shooting')

    % Generate matrices for transcription method
    data=transcriptionMatrix(options.discretization,nx,nu,tau,data);
    
    % p/hp method (LGR)
    if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')

        % Calculate the index transformation for alternative ordering
        if options.reorderLGR
            [ data ] = ReOrderZLGR( data );
        end
        data.map.LGR=LGR;

%         if ~strcmp(options.derivatives,'analytic')
%             [data.FD.vector,data.FD.index,Jac_templete,data.infoForLinkConst]=getPertubationsLGR(sparsity,data.sizes,data);
%         else
            data_temp=data;
            [data.FD.vector,data.FD.index,Jac_templete]=getPertubationsLGR(sparsity,data.sizes,data);
            
       if isfield(problem,'mpflag')
            data_temp.options.derivatives='numeric';
            data.infoForLinkConst=getPertubationsLGR(sparsity_num,data.sizes,data_temp);
       end
            
%         end

        % Constaint Jacobian Structure
        dfz=[Jac_templete.dfxu repmat(sparsity.dfdp,M,1) repmat(ones(1,nt),M*n,1)]; %FD of variables
        idxrowstart=1;
        idxcolstart=1;
        D_structure=zeros(M,M+1); %Structure of the Radal differentiation matrix (multi-segment formulation)
        w=zeros(M,1);
        for i=1:nps
              D_structure(idxrowstart:idxrowstart+npd(i)-1,idxcolstart:idxcolstart+npd(i))=LGR.diff_matrix{npduidx(i)};
              w(idxrowstart:idxrowstart+npd(i)-1)=LGR.weights{npduidx(i)};
              idxrowstart=idxrowstart+npd(i);
              idxcolstart=idxcolstart+npd(i);
        end
        data.map.D_structure=D_structure;
        data.map.w=w;
        
        dfz_noD=dfz;
        dfz=dfz+[kron(speye(n),D_structure) zeros(M*n,M*m+np+nt)];
        dfz_Hrm=[repmat(D_structure,n,n) repmat(D_structure(1:M,1:M),n,m) sparse(n*(M),np+nt) ; repmat(D_structure,m,n) repmat(D_structure(1:M,1:M),m,m) sparse(m*M,np+nt)];
        if ng
            dgz=[Jac_templete.dgxu repmat(sparsity.dgdp,M,1) repmat(ones(1,nt),M*ng,1)];
        else
            dgz=sparse(ng,nt+np+(M+1)*n+M*m);
        end
        
        if nrc
            drcz=[blkdiag(kron(speye(n),D_structure),kron(speye(m),D_structure(:,1:end-1))) zeros(M*(m+n),np) ones(M*(m+n),nt)];
            idx_st_rcl=([rcxl rcul+n]-1)*M+1;
            idx_st_rcu=([rcxu rcuu+n]-1)*M+1;
            idx_st_rce=([rcxe rcue+n]-1)*M+1;

            drcz_rcl=drcz(linspaceMat( idx_st_rcl',idx_st_rcl'+M-1,M ),:);
            drcz_rcu=drcz(linspaceMat( idx_st_rcu',idx_st_rcu'+M-1,M ),:);
            drcz_rce=drcz(linspaceMat( idx_st_rce',idx_st_rce'+M-1,M ),:);
            
            data = transcriptionMatrixRC( data );
            

        else
            drcz_rcl=[];
            drcz_rcu=[];
            drcz_rce=[];
        end
        
        data.map.rcl=rcl;
        data.map.rcu=rcu;
            
        if nb
            dbz=zeros(nb,nz);
            dbz(:,data.FD.index.b(1,:))=1;
        else
            dbz=sparse(nb,nt+np+(M+1)*n+M*m);
        end
        
        data.map.bl=bl;
        data.map.bu=bu;
        data.map.bTol=bTol;

        data.map.spmatsize.jSf=nnz(dfz);
        data.map.spmatsize.jSg=nnz(dgz(gAllidx,:));
        data.map.spmatsize.jSrc=nnz([drcz_rcl;drcz_rcu;drcz_rce]);
        data.map.spmatsize.jSb=nnz(dbz);
            
        jS= [dfz;dgz(gAllidx,:);drcz_rcl;drcz_rcu;drcz_rce;dbz]; %Jacobian structure
        jS_noB=[dfz_noD;dgz(gAllidx,:);zeros(size(drcz_rcl));zeros(size(drcz_rcu));zeros(size(drcz_rce));dbz]; %Jacobian structure excluding the Radau differentiation matrix
        jS_Hrm=[dfz_Hrm;dgz(gAllidx,:);zeros(size(drcz_rcl));zeros(size(drcz_rcu));zeros(size(drcz_rce));dbz]; %Jacobian structure excluding the Radau differentiation matrix
        data.jacStruct=spones(jS);
        data.jS_noB=spones(jS_noB);
        if strcmp(options.transcription,'integral_res_min') || strcmp(options.discretization,'resMinInterpolationForSolution') 
            data.jacStruct_resmin=[data.jacStruct(n*M+1:end,:); sparse(ones(n+ng_eq+1,length(infoNLP.zl)))];
        end
        if isfield(problem,'mpflag') && problem.mpflag
            data.jSidx.org.XUg.row=1:size(jS,1)-nb;
            data.jSidx.org.XUg.col=1:nz-nt-np;
            data.jS_XUg=spones(jS(data.jSidx.org.XUg.row,data.jSidx.org.XUg.col));
            data.jS_XUg_noB=spones(jS_noB(data.jSidx.org.XUg.row,data.jSidx.org.XUg.col));
            if np
                data.jSidx.org.P.row=1:size(jS,1)-nb;
                data.jSidx.org.P.col=nz-nt-np+1:nz-nt;
                data.jS_P=spones(jS(data.jSidx.org.P.row,data.jSidx.org.P.col));
            else
                data.jSidx.org.P.row=[];
                data.jSidx.org.P.col=[];
                data.jS_P=spones(zeros(size(jS,1)-nb,0));
            end
            if nt
                data.jSidx.org.T.row=1:size(jS,1)-nb;
                data.jSidx.org.T.col=nz-nt+1:nz;
                data.jS_T=spones(jS(data.jSidx.org.T.row,data.jSidx.org.T.col));
            else
                data.jSidx.org.T.row=[];
                data.jSidx.org.T.col=[];
                data.jS_T=spones(zeros(size(jS,1)-nb,0));
            end
            if nb
                data.jSidx.org.B_XUg.row=size(jS,1)-nb+1:size(jS,1);
                data.jSidx.org.B_XUg.col=1:nz-nt-np;
                data.jS_B_XUg=spones(jS(data.jSidx.org.B_XUg.row,data.jSidx.org.B_XUg.col));
                if np
                    data.jSidx.org.B_P.row=size(jS,1)-nb+1:size(jS,1);
                    data.jSidx.org.B_P.col=nz-nt-np+1:nz-nt;
                    data.jS_B_P=spones(jS(data.jSidx.org.B_P.row,data.jSidx.org.B_P.col));
                else
                    data.jSidx.org.B_P.row=[];
                    data.jSidx.org.B_P.col=[];
                    data.jS_B_P=spones(zeros(nb,0));
                end
                if nt
                    data.jSidx.org.B_T.row=size(jS,1)-nb+1:size(jS,1);
                    data.jSidx.org.B_T.col=nz-nt+1:nz;
                    data.jS_B_T=spones(jS(data.jSidx.org.B_T.row,data.jSidx.org.B_T.col));
                else
                    data.jSidx.org.B_T.row=[];
                    data.jSidx.org.B_T.col=[];
                    data.jS_B_T=spones(zeros(nb,0));
                end
            else
                data.jSidx.org.B_XUg.row=[];
                data.jSidx.org.B_XUg.col=[];
                data.jS_B_XUg=spones(zeros(0,nz-nt-np));
                data.jSidx.org.B_P.row=[];
                data.jSidx.org.B_P.col=[];
                data.jS_B_P=spones(zeros(nb,0));
                data.jSidx.org.B_T.row=[];
                data.jSidx.org.B_T.col=[];
                data.jS_B_T=spones(zeros(nb,0));
            end
        end
        % Cost (Legrange and mayer) Structure
        dLxu=sparse(M,(M+1)*n+M*m);
        for Li=1:size(sparsity.dLdx,1)
            for xi=1:size(sparsity.dLdx,2)
                dLxu(((Li-1)*M+1):Li*M,((xi-1)*(M+1)+1):(xi*(M+1)-1))=kron(speye(M),sparsity.dLdx(Li,xi));
            end
            for ui=1:size(sparsity.dLdu,2)
                dLxu(((Li-1)*M+1):Li*M,(n*(M+1)+(ui-1)*M+1):(n*(M+1)+ui*M))=kron(speye(M),sparsity.dLdu(Li,ui));
            end
        end
        Lz=[dLxu repmat(sparsity.dLdp,M,1) repmat(ones(1,nt),M,1)];
        
        Ez=sparse(1,unique(data.FD.index.Ey),1,1,nt+np+(M+1)*n+M*m);
        [ib,~,sb]=find(data.FD.index.b);
        data.costStruct.B=sparse(ib,sb,1,nb,nt+np+(M+1)*n+M*m);
        data.costStruct.E=Ez;
        data.costStruct.L=Lz;

        % Legrange Hessian Structure
        data.hessianStruct=spalloc((M+1)*n+M*m,(M+1)*n+M*m,M*(n+m)*(n+m));

        data.map.spmatsize.hSL=nnz(Lz'*Lz);
        data.map.spmatsize.hSE=nnz(Ez'*Ez);
        data.map.spmatsize.hSf=nnz(jS_noB'*jS_noB);
        data.map.spmatsize.hSg=nnz(jS_noB'*jS_noB);

        data.hessianStruct=Lz'*Lz+Ez'*Ez+(jS_noB'*jS_noB);
        if strcmp(options.transcription,'integral_res_min') || strcmp(options.discretization,'resMinInterpolationForSolution')
            data.hessianStruct_resmin=Lz'*Lz+Ez'*Ez+(jS_Hrm'*jS_Hrm);
        elseif strcmp(options.transcription,'direct_collocation_intres_reg')
            data.hessianStruct=Lz'*Lz+Ez'*Ez+(jS_noB'*jS_noB)+(jS_Hrm'*jS_Hrm);
        end
        if options.reorderLGR
             data.hessianStruct=data.hessianStruct(data.reorder.z_idx,data.reorder.z_idx);
        end
        data.hessianStruct=tril(data.hessianStruct);
        if strcmp(options.transcription,'integral_res_min') || strcmp(options.discretization,'resMinInterpolationForSolution')
            data.hessianStruct_resmin=tril(data.hessianStruct_resmin);
        end
        
        hS=data.hessianStruct;
  
        if isfield(problem,'mpflag') && problem.mpflag
            data.hSidx.org.XUXU.row=1:nz-nt-np;
            data.hSidx.org.XUXU.col=1:nz-nt-np;
            data.hS_XUXU=spones(hS(data.hSidx.org.XUXU.row,data.hSidx.org.XUXU.col));
            if np
                data.hSidx.org.PXU.row=nz-nt-np+1:nz-nt;
                data.hSidx.org.PXU.col=1:nz-nt-np;
                data.hS_PXU=spones(hS(data.hSidx.org.PXU.row,data.hSidx.org.PXU.col));
                data.hSidx.org.PP.row=nz-nt-np+1:nz-nt;
                data.hSidx.org.PP.col=nz-nt-np+1:nz-nt;
                data.hS_PP=spones(hS(data.hSidx.org.PP.row,data.hSidx.org.PP.col));
            else
                data.hSidx.org.PXU.row=[];
                data.hSidx.org.PXU.col=[];
                data.hS_PXU=zeros(0,nz-nt-np);
                data.hSidx.org.PP.row=[];
                data.hSidx.org.PP.col=[];
                data.hS_PP=zeros(0,0);
            end
            if nt
                data.hSidx.org.TXU.row=nz-nt+1:nz;
                data.hSidx.org.TXU.col=1:nz-nt-np;
                data.hS_TXU=spones(hS(data.hSidx.org.TXU.row,data.hSidx.org.TXU.col));
                data.hSidx.org.TT.row=nz-nt+1:nz;
                data.hSidx.org.TT.col=nz-nt+1:nz;
                data.hS_TT=spones(hS(data.hSidx.org.TT.row,data.hSidx.org.TT.col));
            else
                data.hSidx.org.TXU.row=[];
                data.hSidx.org.TXU.col=[];
                data.hS_TXU=zeros(0,nz-nt-np);
                data.hSidx.org.TT.row=[];
                data.hSidx.org.TT.col=[];
                data.hS_TT=zeros(0,0);
            end
        end
        
        data.funcs.hessianstructure  = @hessianstructure;
        data.funcs.hessian           = @computeHessian;

    else % h methods
                
%         if ~strcmp(options.derivatives,'analytic')
%             [data.FD.vector,data.FD.index,data.infoForLinkConst]=getPertubations(sparsity,data.sizes,data);
%         else
            data_temp=data;
            [data.FD.vector,data.FD.index]=getPertubations(sparsity,data.sizes,data);
       if isfield(problem,'mpflag')
            data_temp.options.derivatives='numeric';
            data.infoForLinkConst=getPertubations(sparsity_num,data.sizes,data_temp);
       end
%         end
        
        if nrc
            data = transcriptionMatrixRC( data );
        end
        % Exploit structure to generate FD vectors
        %---------------------------------------

            %
            % Generate structure of Jacobian for the constraints for direct
            % collocations methods
            %---------------------------------------
            %

            Fxu=[kron(speye(M/N),sparsity.dfdx),repmat(sparsity.dfdu,M/N,1)];
%             dfz=[rand((M)*n,nt)+1 rand((M)*n,np)+1 kron(speye(N),Fxu)];
            seq=1:(M*n);
            dfz=[repmat(seq',1,nt) repmat(seq',1,np) kron(speye(N),Fxu)];

            Gxu=[kron(speye((M)/N),sparse(sparsity.dgdx)),repmat(sparse(sparsity.dgdu),(M)/N,1)];
            dgz=[repmat(sparsity.dgdt,M,1) repmat(sparsity.dgdp,M,1) kron(speye(N),Gxu)];

            RCxu_1=[repmat(sparse(sparsity.drcdx),(M)/N,ns),repmat(sparse(sparsity.drcdu),(M)/N,ns)];
            RCxu_2=[zeros(size(sparsity.drcdx,1)*(M)/N,m*(ns-1)+n*(ns-1)),repmat(sparse(sparsity.drcdx),(M)/N,1),repmat(sparse(sparsity.drcdu),(M)/N,1)];
            RCxu_1=kron(speye((M-1)/ns),RCxu_1);
            RCxu_2=kron(speye((M-1)/ns),RCxu_2);
            RCxu_1=[RCxu_1 zeros(size(RCxu_1,1),m+n)];
            RCxu_2=[zeros(size(RCxu_2,1),(m+n)) RCxu_2];
            RCxu=RCxu_1+RCxu_2;
            drcz=[repmat(sparsity.drcdt,(M-1)/ns,1) repmat(sparsity.drcdp,(M-1)/ns,1)];
            drcz=[drcz RCxu];
            idx_rcl=logical(repmat([ones(1,vnrcl),zeros(1,vnrcu),zeros(1,vnrce)],1,(M-1)/ns));
            idx_rcu=logical(repmat([zeros(1,vnrcl),ones(1,vnrcu),zeros(1,vnrce)],1,(M-1)/ns));
            idx_rce=logical(repmat([zeros(1,vnrcl),zeros(1,vnrcu),ones(1,vnrce)],1,(M-1)/ns));
            
            A=data.map.A;B=data.map.B;Vx=data.map.Vx;
            data.map.rcl=rcl;
            data.map.rcu=rcu;
            data.map.bl=bl;
            data.map.bu=bu;
            data.map.bTol=bTol;
        
            if N==1
                data.map.spmatsize.jSf=nnz([[sparse(n,nt), sparse(n,np), speye(n), sparse(n,nx+nu-n)]*cx0;A*Vx+B*dfz]);
                data.map.spmatsize.jSg=nnz(dgz(gAllidx,:));
                data.map.spmatsize.jSrc=nnz([drcz(idx_rcl,:);drcz(idx_rcu,:);drcz(idx_rce,:)]);
                data.map.spmatsize.jSb=nnz([sparsity.dbdt0 sparsity.dbdtf sparsity.dbdp sparsity.dbdx0 sparse(nb,(((M)/N-1))*n),...
                sparsity.dbdu0 sparse(nb,nx+nu-2*(n+m)-(((M)/N-1))*n),...
                sparsity.dbdxf]);
            
                jS= [[sparse(n,nt), sparse(n,np), speye(n), sparse(n,nx+nu-n)]*cx0;...
                A*Vx+B*dfz;...
                dgz(gAllidx,:);...
                drcz(idx_rcl,:);drcz(idx_rcu,:);drcz(idx_rce,:);...
                [sparsity.dbdt0 sparsity.dbdtf sparsity.dbdp sparsity.dbdx0 sparse(nb,(((M)/N-1))*n),...
                sparsity.dbdu0 sparse(nb,nx+nu-2*(n+m)-(((M)/N-1))*n),...
                sparsity.dbdxf]];
            else
                data.map.spmatsize.jSf=nnz([[sparse(n,nt), sparse(n,np), speye(n) sparse(n,nx+nu-n)]*cx0;A*Vx+B*dfz]);
                data.map.spmatsize.jSg=nnz(dgz(gAllidx,:));
                data.map.spmatsize.jSrc=nnz([drcz(idx_rcl,:);drcz(idx_rcu,:);drcz(idx_rce,:)]);
                data.map.spmatsize.jSb=nnz([sparsity.dbdt0 sparsity.dbdtf sparsity.dbdp sparsity.dbdx0 sparse(nb,(((M)/N-1))*n),...
                sparsity.dbdu0 sparse(nb,nx+nu-2*(n+m)-(((M)/N-1))*n),...
                sparsity.dbdxf, sparsity.dbduf]);
                
                jS=[[sparse(n,nt), sparse(n,np), speye(n) sparse(n,nx+nu-n)]*cx0;...
                A*Vx+B*dfz;
                dgz(gAllidx,:);...
                drcz(idx_rcl,:);drcz(idx_rcu,:);drcz(idx_rce,:);...
                [sparsity.dbdt0 sparsity.dbdtf sparsity.dbdp sparsity.dbdx0 sparse(nb,(((M)/N-1))*n),...
                sparsity.dbdu0 sparse(nb,nx+nu-2*(n+m)-(((M)/N-1))*n),...
                sparsity.dbdxf, sparsity.dbduf]]; 
            end
            
           data.jacStruct=spones(jS);
           if strcmp(options.transcription,'integral_res_min') || strcmp(options.discretization,'resMinInterpolationForSolution')
                data.jacStruct_resmin=[data.jacStruct(n*M+1:end,:); sparse(ones(n+ng_eq+1,length(infoNLP.zl)))];
           end
           
           if isfield(problem,'mpflag') && problem.mpflag
                data.jSidx.org.XUg.row=1:size(jS,1)-nb;
                data.jSidx.org.XUg.col=1+nt+np:nz;
                data.jS_XUg=spones(jS(data.jSidx.org.XUg.row,data.jSidx.org.XUg.col));
                if np
                    data.jSidx.org.P.row=1:size(jS,1)-nb;
                    data.jSidx.org.P.col=nt+1:nt+np;
                    data.jS_P=spones(jS(data.jSidx.org.P.row,data.jSidx.org.P.col));
                else
                    data.jSidx.org.P.row=[];
                    data.jSidx.org.P.col=[];
                    data.jS_P=spones(zeros(size(jS,1)-nb,0));
                end
                if nt
                    data.jSidx.org.T.row=1:size(jS,1)-nb;
                    data.jSidx.org.T.col=1:nt;
                    data.jS_T=spones(jS(data.jSidx.org.T.row,data.jSidx.org.T.col));
                else
                    data.jSidx.org.T.row=[];
                    data.jSidx.org.T.col=[];
                    data.jS_T=spones(zeros(size(jS,1)-nb,0));
                end
                if nb
                    data.jSidx.org.B_XUg.row=size(jS,1)-nb+1:size(jS,1);
                    data.jSidx.org.B_XUg.col=1+nt+np:nz;
                    data.jS_B_XUg=spones(jS(data.jSidx.org.B_XUg.row,data.jSidx.org.B_XUg.col));
                    if np
                        data.jSidx.org.B_P.row=size(jS,1)-nb+1:size(jS,1);
                        data.jSidx.org.B_P.col=nt+1:nt+np;
                        data.jS_B_P=spones(jS(data.jSidx.org.B_P.row,data.jSidx.org.B_P.col));
                    else
                        data.jSidx.org.B_P.row=[];
                        data.jSidx.org.B_P.col=[];
                        data.jS_B_P=spones(zeros(nb,0));
                    end
                    if nt
                        data.jSidx.org.B_T.row=size(jS,1)-nb+1:size(jS,1);
                        data.jSidx.org.B_T.col=1:nt;
                        data.jS_B_T=spones(jS(data.jSidx.org.B_T.row,data.jSidx.org.B_T.col));
                    else
                        data.jSidx.org.B_T.row=[];
                        data.jSidx.org.B_T.col=[];
                        data.jS_B_T=spones(zeros(nb,0));
                    end
                else
                    data.jSidx.org.B_XUg.row=[];
                    data.jSidx.org.B_XUg.col=[];
                    data.jS_B_XUg=spones(zeros(0,nz-nt-np));
                    data.jSidx.org.B_P.row=[];
                    data.jSidx.org.B_P.col=[];
                    data.jS_B_P=spones(zeros(nb,0));
                    data.jSidx.org.B_T.row=[];
                    data.jSidx.org.B_T.col=[];
                    data.jS_B_T=spones(zeros(nb,0));
                end
           end
          
           
            %
            % Generate structure of Hessian for the cost function
            %---------------------------------------

            Lxu=kron(speye(N),[kron(speye((M)/N),sparsity.dLdx),repmat(sparsity.dLdu,(M)/N,1)]);
            Lz=[ones(M,nt) repmat(sparsity.dLdp,M,1) Lxu];

            Ez=sparse(1,data.FD.index.Ey,1,1,nt+np+n*M+m*N);
            [ib,~,sb]=find(data.FD.index.b);
            data.costStruct.B=sparse(ib,sb,1,nb,nt+np+n*M+m*N);
            data.costStruct.E=Ez;
            data.costStruct.L=Lz;

%             data.hessianStruct=spalloc(M*(n+m),M*(n+m),M*(n+m)*(n+m));
            
            data.map.spmatsize.hSL=nnz(Lz'*Lz);
            data.map.spmatsize.hSE=nnz(Ez'*Ez);
            data.map.spmatsize.hSf=nnz(data.jacStruct'*data.jacStruct);
            data.map.spmatsize.hSg=nnz(data.jacStruct'*data.jacStruct);
            
            data.hessianStruct=spones(tril(Lz'*Lz+Ez'*Ez+(data.jacStruct'*data.jacStruct)));
            data.hessianStruct_resmin=data.hessianStruct;
            data.funcs.hessianstructure  = @hessianstructure;
            data.funcs.hessian           = @computeHessian;
            
            hS=data.hessianStruct;
            
            if isfield(problem,'mpflag') && problem.mpflag
                data.hSidx.org.XUXU.row=1+nt+np:nz;
                data.hSidx.org.XUXU.col=1+nt+np:nz;
                data.hS_XUXU=spones(hS(data.hSidx.org.XUXU.row,data.hSidx.org.XUXU.col));
                if np
                    data.hSidx.org.PXU.row=1+nt+np:nz;
                    data.hSidx.org.PXU.col=1+nt:nt+np;
                    data.hS_PXU=spones(hS(data.hSidx.org.PXU.row,data.hSidx.org.PXU.col));
                    data.hSidx.org.PP.row=1+nt:nt+np;
                    data.hSidx.org.PP.col=1+nt:nt+np;
                    data.hS_PP=spones(hS(data.hSidx.org.PP.row,data.hSidx.org.PP.col));
                else
                    data.hSidx.org.PXU.row=[];
                    data.hSidx.org.PXU.col=[];
                    data.hS_PXU=zeros(0,nz-nt-np);
                    data.hSidx.org.PP.row=[];
                    data.hSidx.org.PP.col=[];
                    data.hS_PP=zeros(0,0);
                end
                if nt
                    data.hSidx.org.TXU.row=1+nt+np:nz;
                    data.hSidx.org.TXU.col=1:nt;
                    data.hS_TXU=spones(hS(data.hSidx.org.TXU.row,data.hSidx.org.TXU.col));
                    data.hSidx.org.TT.row=1:nt;
                    data.hSidx.org.TT.col=1:nt;
                    data.hS_TT=spones(hS(data.hSidx.org.TT.row,data.hSidx.org.TT.col));
                else
                    data.hSidx.org.TXU.row=[];
                    data.hSidx.org.TXU.col=[];
                    data.hS_TXU=zeros(0,nz-nt-np);
                    data.hSidx.org.TT.row=[];
                    data.hSidx.org.TT.col=[];
                    data.hS_TT=zeros(0,0);
                end
            end
    end

    % Format sparsity pattern for Worhp NLP Solver
    if strcmp(options.NLPsolver,'worhp')
        [data.jacStruct_GRow, data.jacStruct_GCol]=find(data.jacStruct);
        [jacStruct_HMRow, jacStruct_HMCol]=find(tril(data.hessianStruct,-1));
        [jacStruct_HMRowD, jacStruct_HMColD]=find(diag(diag(data.hessianStruct)));
        data.jacStruct_HMRow=[jacStruct_HMRow;jacStruct_HMRowD];
        data.jacStruct_HMCol=[jacStruct_HMCol;jacStruct_HMColD];
    end

else   % if options.transcription,'multiple_shooting'
 
   % Create matrix to extract x vector from z
   Vx=[zeros(nx,nt+np) [kron(speye(N),[speye(n) zeros(n,m)]); zeros(n,(n+m)*N)] [zeros((M-1)*n,n);speye(n)]];
   xV=Vx\speye(n*M);
     
   % Create matrix to extract u vector from z
   Vu=[zeros(nu,nt+np) kron(speye(N),[zeros(m,n) speye(m)])  zeros(nu,n)]; 
   uV=Vu\speye(m*N);
   
   data.map.Vu=Vu; 
   data.map.Vx=Vx;     
   data.map.xV=xV;
   data.map.uV=uV;

   % Exploit structure to generate FD vectors 
   %---------------------------------------

   [data.FD.vector,data.FD.index]=getPertubations_multipleShooting(sparsity,data.sizes,data);  
   Ez=sparse(1,data.FD.index.Ey,1,1,nt+np+n*M+m*N);
   data.costStruct.E=Ez;
   [ib,jb,sb]=find(data.FD.index.b);
   data.costStruct.B=sparse(ib,sb,1,nb,nt+np+n*M+m*N);
    
   
    %
    % Generate structure of Jacobian for the constraints for the
    % multiple-shooting
    %---------------------------------------
    %
         
   dcdz=[zeros(n,nt), zeros(n,np), eye(n), zeros(n,nx+nu-n)]*cx0;     % Dynamic constraints
   dgdz=[];   % Path Constraints
   dbdz=[]; 
   nm=n+m; 
    
       
   if ng
      dgdz=[dgdz; sparsity.dgdt,  sparsity.dgdp, sparsity.dgdx, sparsity.dgdu, zeros(ng,nm*N-m)];
   end


   for i=0:M-3 % Compute costs,constraints and sensitivities M-1 times
     k=i*nm+1;  
        dcdz=[dcdz;zeros(n,nt), sparsity.dfdp, zeros(n,k-1), ones(n,nm), -eye(n), zeros(n,(nm)*(M-i-2))];
      if ng
        dgdz=[dgdz; sparsity.dgdt,  sparsity.dgdp,  zeros(ng,k-1+nm), sparsity.dgdx, sparsity.dgdu, zeros(ng,nm*(N-i-1)-m)];
      end
   end
   dcdz=[dcdz;zeros(n,nt), sparsity.dfdp, zeros(n,(M-2)*nm), ones(n,nm), -eye(n)];



   if nb
     dbdz=[sparsity.dbdtf, sparsity.dbdp, sparsity.dbdx0, sparsity.dbdu0, zeros(nb,nx+nu-2*(n+m)), sparsity.dbduf, sparsity.dbdxf];
     bzz=spones(dbdz'*dbdz);
    else 
     bzz=0;
    end

    jS=[dcdz;dgdz;dbdz];
    data.jacStruct=spones(jS);
   
end

%Offline generation of adigator files
 if strcmp(options.derivatives,'adigator') 
     if adigatorGen==1
        if ~strcmp(options.transcription,'integral_res_min')
            data=genAdigator4ICLOCS( options, data, n, m, np, nt, M );
        end
     else
        currentFolder = pwd;
        cd(options.adigatorPath);
        startupadigator;
        cd(currentFolder);
        if ~strcmp(options.transcription,'integral_res_min')
            data=genAdigator4ICLOCS( options, data, n, m, np, nt, M );
        end
        adigatorGen=1;
     end
 end



% Format initial guess/reference for NLP
%---------------------------------------
infoNLP.z0=zeros(nz,1);
data.data.map_w=data.map.w;

if isempty(guess.states)
   
  x0_lg=x0_l;x0_lg(x0_l==-inf)=-100;
  x0_ug=x0_u;x0_ug(x0_u==inf)=100;
  xf_lg=xf_l;xf_lg(xf_l==-inf)=-100;
  xf_ug=xf_u;xf_ug(xf_u==inf)=100;
  guess.states=[(x0_lg+(x0_ug-x0_lg).*rand(n,1))';(xf_lg+(xf_ug-xf_lg).*rand(n,1))'];    % Prune bounds for initial conditions
end

if isfield(guess,'inputs') && isempty(guess.inputs)
   ul_g=ul;ul_g(ul==-inf)=-100;
   uu_g=uu;uu_g(uu==inf)=100;
   guess.inputs=[(ul_g+(uu_g-ul_g).*rand(m,1))';(ul_g+(uu_g-ul_g).*rand(m,1))'];
end    
 

[ x_guess, u_guess, data ] = getGuessSolution( options, guess, problem, data );

if strcmp(options.transcription,'multiple_shooting')
    
    if size(u_guess,1)/m==size(x_guess,1)/n
        u_guess=u_guess(1:end-m,:)';u_guess=u_guess(:);
    else
        u_guess=u_guess';u_guess=u_guess(:);
    end
    x_guess=x_guess';x_guess=x_guess(:);
    infoNLP.z0=data.map.xV*x_guess+data.map.uV*u_guess;
else
    infoNLP.z0=data.map.xV*x_guess+data.map.uV*u_guess;
end

if isfield(problem.states,'x0') && ~isempty(problem.states.x0)
    if strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR')
        if data.options.reorderLGR
            infoNLP.z0(1:n)=data.x0t;
        else
            infoNLP.z0(1:M+1:M*n)=data.x0t;
        end
    else
        infoNLP.z0(nt+np+1:n+nt+np)=data.x0t;
    end
end

if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
        if np;infoNLP.z0(n*(M+1)+M*m+1:n*(M+1)+M*m+np)=guess.parameters(:);end
        if nt==1
            infoNLP.z0(nz)=guess.tf;
        elseif nt>=2
            infoNLP.z0(nz-nt+1:nz)=linspace(guess.t0,guess.tf,nt);
        end
else
    if np;infoNLP.z0(nt+1:nt+np)=guess.parameters(:);end
    if nt==1
        infoNLP.z0(1)=guess.tf;
    elseif nt==2
        infoNLP.z0(1)=guess.t0;
        infoNLP.z0(2)=guess.tf;
    end
end


  %  Formatting options of the ipopt solver
  data.funcs.objective         = @costFunction;
  data.funcs.gradient          = @costGradient;
  data.funcs.constraints       = @constraintFunction;
  data.funcs.jacobian          = @constraintJacobian;
  data.funcs.jacobianstructure = @jacobianstructure;
  
  if isfield(problem,'callback') && ~isempty(problem.callback)
      data.funcs.iterfunc=problem.callback;
  end
  
  % Set function types 
  [ data ] = setFunctionTypes( problem, data );
  
  %Alternative ordering of LGR method 
  if  ((strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))) && options.reorderLGR
      infoNLP.zl=infoNLP.zl(data.reorder.z_idx);
      infoNLP.zu=infoNLP.zu(data.reorder.z_idx);
      infoNLP.z0=infoNLP.z0(data.reorder.z_idx);
      infoNLP.cl=infoNLP.cl(data.reorder.vert_idx);                    % Lower bounds on constraints.
      infoNLP.cu=infoNLP.cu(data.reorder.vert_idx);                      % Upper bounds on constraints
      data.jacStruct=data.jacStruct(data.reorder.vert_idx,data.reorder.z_idx);
  end
  
  % Separate Ceq and Cneq
  
    infoNLP.nnod=n*M;
    flag_eq=(infoNLP.cu(infoNLP.nnod+1:end)-infoNLP.cl(infoNLP.nnod+1:end)==0);
    flag_inf_ub=isinf(infoNLP.cu(infoNLP.nnod+1:end));
    flag_inf_lb=isinf(infoNLP.cl(infoNLP.nnod+1:end));
    infoNLP.ind_eqODE=1:infoNLP.nnod;
    infoNLP.ind_eq=find(flag_eq); 
    infoNLP.ind_ineq_ub=find(~flag_eq & ~flag_inf_ub);
    infoNLP.ind_ineq_lb=find(~flag_eq & ~flag_inf_lb);

  
  data.infoNLP=infoNLP;

  if strcmp(options.transcription,'integral_res_min')
    if  ((strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR')))
        t_list=(data.tau_segment+1)/2;
    else
        t_list=[0;data.tau_inc(2:ns:end)]/ns;
    end
    data.data.resmin=1;
    data = transcribeResMin( t_list,options,data );
    
    
    data.funcs=data.dataNLP.funcs;
    data.funcs.jacobianstructure = @jacobianstructure_resmin;
%     if  ((strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR')))
%         data.funcs.hessianstructure  = @(data) sparse(tril(ones(length(infoNLP.z0),length(infoNLP.z0))));
%     else
        data.funcs.hessianstructure  = @hessianstructure_resmin;
%     end
    
    data.options.transcription= data.dataNLP.options.transcription;
    
    if strcmp(options.derivatives,'adigator')
        genAdigator4ICLOCS_resmin( options, data, n, m, np, nt, M ,ng_eq );
    end
 
  else
    if  ((strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR')))
        t_list=(data.tau_segment'+1)/2;
    else
        t_list=[0;data.tau_inc(2:ns:end)]/ns;
    end
    if strcmp(options.transcription,'direct_collocation') && strcmp(options.discretization,'resMinInterpolationForSolution')
        data.data.resmin=1;
    else
        data.data.resmin=0;
    end
    
    if strcmp(options.transcription,'direct_collocation_intres_reg')
        data.resmin = transcribeResMin( t_list,options,data );
        if strcmp(options.derivatives,'adigator')
            genAdigator4ICLOCS_resmin( options, data.resmin, n, m, np, nt, M ,ng_eq );
        end
    elseif isfield(options.print,'residual_error') && options.print.residual_error
        data.resmin = transcribeResErrorAnalysis( t_list,options,data );
    end
    
  end

  
  

 
% Check syntax of user defined dynamics and path constraints.
data=checkDynamics( infoNLP.z0,data );




end

