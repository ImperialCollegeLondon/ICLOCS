function [infoNLP,data,options]=transcribeOCP(problem,guess,options)
%TRANSCRIBEOCP - Process information from 'problem', 'guess' and 'options' for NLP solver
%Specifically:
%Error checking of function definitions and bounds
%Define bounds for NLP variable + continuity, path and boundary constraints
%Format matrices for direct transcription method(if required)
%Generate initial guess for optimization
%Generate structure of the jacobian of the constraints(if required)
%Construct optimal finite-difference pertubation sets(if required)
%
% Syntax:  [infoNLP,data]=transcribeOCP(problem,guess,options)
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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk


% Discretization Error Tol
options.discErrorTol=problem.states.xErrorTol;

% Constraint Violation Tol
options.constraintErrorTol=problem.constraintErrorTol;

% Automatic discretization 
if (strcmp(options.transcription,'AutoDirect')) && ~isfield(options,'AutoDirect')
    options.AutoDirect=1;
    options.transcription='hermite';
end

if strcmp(options.transcription,'hpLGR') && options.adaptseg==1
    problem.constraints.bu=[problem.constraints.bl,ones(1,options.nsegment)*options.maxtimeinterval];
    problem.constraints.bl=[problem.constraints.bl,ones(1,options.nsegment)*options.mintimeinterval];
end

% Define (and assign) some parameters
%---------------------------------------
n=length(problem.states.x0l);              % Number of states
m=length(problem.inputs.ul);               % Number of inputs
np=length(problem.parameters.pl);          % Number of free parameters
ng=length(problem.constraints.gl);         % Number of path constraints
nb=length(problem.constraints.bl);         % Number of boundary constraints
N=problem.inputs.N;                        % Number of control actions
data.data=problem.data;
data.data.transcription=options.transcription;


% Scale problem variables
if options.scaling
    [ problem,guess ] = scale_problem( problem,guess );
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
if (strcmp(options.transcription,'globalLGR')) || (strcmp(options.transcription,'hpLGR'))
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
else
    if length(options.tau)~=1
        M=length(options.tau)+1;
    else    
        M=options.nodes;  
    end
end


% Multiple-shooting: Number of integration nodes in the interval t=[t0,tf] is N+1
if strcmp(options.transcription,'multiple_shooting')
   if N==0
       N=M-1;
   else
       M=N+1;
   end % Default N=M-1.
elseif strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')
    if N==0;N=M-1;end % Default N=M-1.
else
    if N==0;N=M;end % Default N=M-1.
end


% The final time for discrete time systems is imposed equal to 1 in order 
% to use a unform formulation of the optimization problem
if (strcmp(options.transcription,'discrete'))
   problem.time.tf_min=1;
   problem.time.tf_max=1;
end 
 

% Get bounds for final time and check if time is free or fixed
% nt=0 when the final time is not a variable otherwise nt=1.
% For LGR, when t0 and tf both variables: nt=2
%          when adaptive method, nt equal to the actual number of time
%          varibles (>2)
tfl=problem.time.tf_min; tfu=problem.time.tf_max;
if tfl==tfu&&~strcmp(options.transcription,'globalLGR')&&~strcmp(options.transcription,'hpLGR')
    data.tf=tfl; nt=0;tfl=[];tfu=[];
else
    if strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')
        if options.adaptseg==1
            nt=nps+1;
        else
            nt=2;
        end
    else
        nt=1; 
    end
end

%Offline generation of adigator files
 if strcmp(options.derivatives,'adigator')
    genAdigator4ICLOCS( options, data, n, m, np );
 end

% Set the initial time and other parameters to adjust the temporal
% scale
data.t0=problem.time.t0;
if (strcmp(options.transcription,'discrete'))
     data.k0=problem.time.t0;
     problem.time.t0=0;
     data.Nm=N;
elseif (strcmp(options.transcription,'globalLGR')) || (strcmp(options.transcription,'hpLGR'))
     data.Nm=1;
else
     data.k0=problem.time.t0;
     data.Nm=1;
end

% Number of algebraic variable rate constraints
if (strcmp(options.transcription,'hermite'))
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
        
        vnrcl=(sum(~isinf(problem.states.xrl))-sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru))))*3+sum(~isinf(problem.inputs.url))*2;
        vnrcu=(sum(~isinf(problem.states.xru))-sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru))))*3+sum(~isinf(problem.inputs.uru))*2;
        vnrce=sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru)))*2+sum(problem.inputs.url(~isinf(problem.inputs.url))==problem.inputs.uru(~isinf(problem.inputs.uru)))*2;
        nrcl=(sum(~isinf(problem.states.xrl))-sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru))))*3*(M-1)+sum(~isinf(problem.inputs.url))*2*(N-1);
        nrcu=(sum(~isinf(problem.states.xru))-sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru))))*3*(M-1)+sum(~isinf(problem.inputs.uru))*2*(N-1);
        nrce=sum(problem.states.xrl(~isinf(problem.states.xrl))==problem.states.xru(~isinf(problem.states.xru)))*2*(M-1)+sum(problem.inputs.url(~isinf(problem.inputs.url))==problem.inputs.uru(~isinf(problem.inputs.uru)))*2*(N-1);
    else
        vnrcl=0;vnrcu=0;vnrce=0;nrcl=0;nrcu=0;nrce=0;
    end
    M=2*M-1;N=N*2-1;ns=2;
elseif strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')
    ns=1;
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
    ns=1;
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

% Number of path constraints
% Save the corrsponding variables sizes
if (strcmp(options.transcription,'globalLGR')) || (strcmp(options.transcription,'hpLGR'))
    data.sizes={nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce};
else    
    data.sizes={nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce};
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
if isempty(problem.states.x0)
  cx0=0;
  data.x0=x0l;
  data.x0t=(x0_l+(x0_u-x0_l).*rand(n,1));
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
else 
    gl=[];
    gu=[]; 
end

% Get bounds for the constraint functions
if nrc
    rcl=[zeros(1,nrcl),-inf*ones(1,nrcu),zeros(nrce)];
    rcu=[inf*ones(1,nrcl),zeros(1,nrcu),zeros(nrce)];
else 
    rcl=[];
    rcu=[];
end

if nb
    bl=problem.constraints.bl(:);
    bu=problem.constraints.bu(:);
else 
    bl=[];
    bu=[]; 
end


% Store some matrices in data structure
data.options=options;
data.functions=problem.functions;
data.functions_unscaled=problem.functions_unscaled;
if options.scaling
    data.data.Xscale=problem.states.scales;
    data.data.Uscale=problem.inputs.scales;
    data.data.Tscale=problem.time.scales;
    data.data.Tshift=problem.time.shifts;
    data.data.Xshift=problem.states.shifts;
    data.data.Ushift=problem.inputs.shifts;
    if np
        data.data.Pscale=problem.parameters.scales;
        data.data.Pshift=problem.parameters.shifts;
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
data.references.xr=repmat(problem.setpoints.states,M,1);   
data.references.ur=repmat(problem.setpoints.inputs,M,1);   


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
        
    case{'discrete','euler','trapezoidal','hermite'}
        if mod(M,N)||(N>M)              % Check for errors
            error('# integration steps is not divisible or less than N');
        end
%         nk=M-1;                      % Number of integration steps
        nx=M*n;                      % Number of unknown states  
        nu=N*m;                      % Number of unknown controls
	    nz=nt+np+nx+nu;              % Length of the optimization variable

        xpl=repmat(xl(:),M/N,1);xpu=repmat(xu(:),M/N,1);
       
        infoNLP.zl=[tfl; pl;repmat([xpl(:);ul(:)],N-1,1);xf_l(:);ul(:)];
        infoNLP.zu=[tfu; pu;repmat([xpu(:);uu(:)],N-1,1);xf_u(:);uu(:)];
        infoNLP.zl(nt+np+1:nt+np+n)=x0_l;    % Prune bounds for initial conditions
        infoNLP.zu(nt+np+1:nt+np+n)=x0_u;
        infoNLP.zl(nt+np+n+1:nt+np+n+m)=u0_l;    % Prune bounds for initial conditions
        infoNLP.zu(nt+np+n+1:nt+np+n+m)=u0_u;
        
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
            tl_mat=[data.t0 cumsum(ones(1,nt-1)*options.mintimeinterval)];
            tu_mat=[data.t0 cumsum(ones(1,nt-1)*options.maxtimeinterval)];
        else
            tl_mat=[data.t0 tfl];
            tu_mat=[data.t0 tfu];
        end       
        infoNLP.zl=[xl_mat(:);ul_mat(:);pl;tl_mat']; %lower bounds of z
        infoNLP.zu=[xu_mat(:);uu_mat(:);pu;tu_mat']; %upper bounds of z

    otherwise;disp('Unknown method. Check spelling');
end




% Define bounds for constraint functions
%---------------------------------------
% Constraints: c = [c0... cF g(0)....g(ng-1) bo, ,bo(nb-1)]' 
%              ck -> defects  for k=0, ..., F=M-1 
%              gk -> general path constraints g(x(k),u(k),p,k) for
%              k=0,...,ng-1
%              bo(k) -> nb boundary conditions b(x(0),x(f),u(0),u(f),p,t)
if strcmp(options.transcription,'multiple_shooting')
infoNLP.cl=[kron(ones(M,1),zeros(n,1));kron(ones(M-1,1),gl(:));bl(:)];
infoNLP.cu=[kron(ones(M,1),zeros(n,1));kron(ones(M-1,1),gu(:));bu(:)];   
else
infoNLP.cl=[kron(ones(M,1),-eps*ones(n,1));kron(ones(M,1),gl(:));rcl(:);bl(:)];
infoNLP.cu=[kron(ones(M,1),eps*ones(n,1));kron(ones(M,1),gu(:));rcu(:);bu(:)];
end

% Extract sparsity structures
%---------------------------------------
if ~strcmp(options.derivatives,'analytic')
    sparsity=getStructure(data.sizes,options.transcription);
else
    % The structure of the derivatives is determined only considering a
    % a fixed structure
    sparsity=getStructureA(data);
end
data.sparsity=sparsity; 


% Generation of mesh (time dimension)
if strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')
    if isfield(options,'tau_segment')
        tau_segment=options.tau_segment; %user specified mesh segments
    else
        tau_segment=-1:2/nps:1; %uniform mesh segments
    end
    %Calculate all the LGR points, weights and differentiation matrices
    LGR=generateLGR_All(npdu); 
    tau_local_seg=zeros(M,1); %mesh in local (segmental) time frame
    idxtemp=[0,cumsum(npd)];
    for i=1:length(npduidx)
        tau_local_seg(idxtemp(i)+1:idxtemp(i+1),:)=LGR.points{npduidx(i)};
    end
    tau_inc=zeros(M,1); %mesh in incremental (global) time frame
    tau_seg=cell(length(tau_segment)-1,1);
    for i=1:length(tau_segment)-1
        tau_seg{i} = ((tau_segment(i+1)+tau_segment(i)) + (tau_segment(i+1)-tau_segment(i)) * [LGR.points{npduidx(i)}; 1]) / (1-(-1));
        tau_inc(idxtemp(i)+1:idxtemp(i+1),:)=tau_seg{i}(1:end-1);
    end
    tau=diff(tau_inc); %mesh in difference formulation
    data.tau_inc=tau_inc;
    data.tau_local_seg=tau_local_seg;
    data.tau_segment=tau_segment;
else
    if options.tau==0
        tau=ns*ones(M-1,1)/(M-1); 
    else
       texst=ones(ns,1);
       tau=kron(options.tau,texst);
    end
    data.tau_inc=cumsum(tau);
    if abs(sum(tau)-ns)>sqrt(eps);error('Time vector (tau) should sum to 1');end
end
data.tau=tau;


% Format direct transcription matrices
%---------------------------------------
% Format matrices for continuity constraints: c(z)=[A.Vx.z+B.F(z)]
% and generate the quadrature vector w(tau)

if ~strcmp(options.transcription,'multiple_shooting')

    % Generate matrices for transcription method
    data=transcriptionMatrix(options.transcription,nx,nu,tau,data);
    
    % p/hp method (LGR)
    if strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')

        % Calculate the index transformation for alternative ordering
        if options.reorderLGR
            [ data ] = ReOrderZLGR( data );
        end
        data.map.LGR=LGR;
        
        [data.FD.vector,data.FD.index,Jac_templete]=getPertubationsLGR(sparsity,data.sizes,data);

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
        
        if nb
            dbz=zeros(nb,nz);
            dbz(:,data.FD.index.b(1,:))=1;
        else
            dbz=sparse(nb,nt+np+(M+1)*n+M*m);
        end

        jS= [dfz;dgz;drcz_rcl;drcz_rcu;drcz_rce;dbz]; %Jacobian structure
        jS_noB=[dfz_noD;dgz;dbz]; %Jacobian structure excluding the Radau differentiation matrix
        data.jacStruct=spones(jS);
        data.jS_noB=spones(jS_noB);

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
        data.hessianStruct=Lz'*Lz+Ez'*Ez+(jS_noB'*jS_noB);
        if options.reorderLGR
             data.hessianStruct=data.hessianStruct(data.reorder.z_idx,data.reorder.z_idx);
        end
        data.hessianStruct=tril(data.hessianStruct);
        data.funcs.hessianstructure  = @hessianstructure;
        data.funcs.hessian           = @computeHessian;

    else % h methods
        
        [data.FD.vector,data.FD.index]=getPertubations(sparsity,data.sizes,data);

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
            dfz=[rand((M)*n,nt)+1 repmat(sparsity.dfdp,(M),1) kron(speye(N),Fxu)];

            Gxu=[kron(speye((M)/N),sparse(sparsity.dgdx)),repmat(sparse(sparsity.dgdu),(M)/N,1)];
            dgz=[repmat(sparsity.dgdt,M,1) repmat(sparsity.dgdp,M,1) kron(speye(N),Gxu)];

            RCxu_1=[repmat(sparse(sparsity.drcdx),(M)/N,ns),repmat(sparse(sparsity.drcdu),(M)/N,ns)];
            RCxu_2=[zeros(size(sparsity.drcdx,1)*(M)/N,m*(ns-1)+n*(ns-1)),repmat(sparse(sparsity.drcdx),(M)/N,1),repmat(sparse(sparsity.drcdu),(M)/N,1)];
            RCxu_1=kron(speye((M-1)/ns),RCxu_1);
            RCxu_2=kron(speye((M-1)/ns),RCxu_2);
            RCxu_1=[RCxu_1 zeros(size(RCxu_1,1),m+n)];
            RCxu_2=[zeros(size(RCxu_2,1),(m+n)) RCxu_2];
            RCxu=RCxu_1(1:(M-1)/ns:size(RCxu_1,1),:)+RCxu_2(1:(M-1)/ns:size(RCxu_2,1),:);
            drcz=[repmat(sparsity.drcdt,(M-1)/ns,1) repmat(sparsity.drcdp,(M-1)/ns,1)];
            drcz=drcz(1:(M-1)/ns:size(drcz,1),:);
            drcz=[drcz RCxu];
            idx_rcl=logical(repmat([ones(1,vnrcl),zeros(1,vnrcu),zeros(1,vnrce)],1,(M-1)/ns));
            idx_rcu=logical(repmat([zeros(1,vnrcl),ones(1,vnrcu),zeros(1,vnrce)],1,(M-1)/ns));
            idx_rce=logical(repmat([zeros(1,vnrcl),zeros(1,vnrcu),ones(1,vnrce)],1,(M-1)/ns));
            
            A=data.map.A;B=data.map.B;Vx=data.map.Vx;

            if N==1
                jS= [[zeros(n,nt), zeros(n,np), eye(n), zeros(n,nx+nu-n)]*cx0;...
                A*Vx+B*dfz;...
                dgz;...
                drcz(idx_rcl,:);drcz(idx_rcu,:);drcz(idx_rce,:);...
                [sparsity.dbdtf sparsity.dbdp sparsity.dbdx0 zeros(nb,(((M)/N-1))*n),...
                sparsity.dbdu0 zeros(nb,nx+nu-2*(n+m)-(((M)/N-1))*n),...
                sparsity.dbdxf]];
            else
                jS=[[zeros(n,nt), zeros(n,np), eye(n) zeros(n,nx+nu-n)]*cx0;...
                A*Vx+B*dfz;
                dgz;...
                drcz(idx_rcl,:);drcz(idx_rcu,:);drcz(idx_rce,:);...
                [sparsity.dbdtf sparsity.dbdp sparsity.dbdx0 zeros(nb,(((M)/N-1))*n),...
                sparsity.dbdu0 zeros(nb,nx+nu-2*(n+m)-(((M)/N-1))*n) sparsity.dbduf,...
                sparsity.dbdxf]]; 
            end
            
           data.jacStruct=spones(jS);
           
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

            data.hessianStruct=spalloc(M*(n+m),M*(n+m),M*(n+m)*(n+m));
            data.hessianStruct=tril(Lz'*Lz+Ez'*Ez+(data.jacStruct'*data.jacStruct));
            data.funcs.hessianstructure  = @hessianstructure;
            data.funcs.hessian           = @computeHessian;
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

% Format initial guess/reference for NLP
%---------------------------------------
infoNLP.z0=zeros(nz,1);
if strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')
    x_guess=zeros(n,M+1);u_guess=zeros(m,M);
else
    x_guess=zeros(M,n);u_guess=zeros(N,m);
end

if isempty(guess.states)
   
  x0_lg=x0_l;x0_lg(x0_l==-inf)=-100;
  x0_ug=x0_u;x0_ug(x0_u==inf)=100;
  xf_lg=xf_l;xf_lg(xf_l==-inf)=-100;
  xf_ug=xf_u;xf_ug(xf_u==inf)=100;
  guess.states=[(x0_lg+(x0_ug-x0_lg).*rand(n,1))';(xf_lg+(xf_ug-xf_lg).*rand(n,1))'];    % Prune bounds for initial conditions
end

if isempty(guess.inputs)
   ul_g=ul;ul_g(ul==-inf)=-100;
   uu_g=uu;uu_g(uu==inf)=100;
   guess.inputs=[(ul_g+(uu_g-ul_g).*rand(m,1))';(ul_g+(uu_g-ul_g).*rand(m,1))'];
end    
 
% warm starting from guess
if strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')
    if (strcmp(options.start,'Cold'))
        
        if isfield(guess,'time')
            Tx=linspace(0,guess.time(end), M+1);
            Tu=linspace(0,guess.time(end), M);
            x_guess=interp1(guess.time, guess.states,Tx,'linear');
            x_guess=x_guess(:);
            u_guess=interp1(guess.time, guess.inputs,Tu,'linear');
            u_guess=u_guess(:);
        else
            for i=1:n
              x_guess(i,:)=linspace(guess.states(1,i),guess.states(2,i),M+1);
            end
            x_guess=reshape(x_guess',(M+1)*n,1);
            for i=1:m
             if M>1  
               u_guess(i,:)=linspace(guess.inputs(1,i),guess.inputs(2,i),M);
              else
               u_guess(i,:)=guess.inputs(1,i);
              end
            end
            u_guess=reshape(u_guess',M*m,1);
        end
        data.multipliers=[];
    elseif (strcmp(options.start,'Warm'))
        T=(guess.tf-data.t0)/2*data.tau_inc+(guess.tf+data.t0)/2;
        if size(guess.time,1)==size(guess.states,1)
            x_guess=interp1(guess.time, guess.states,[T;guess.tf],'pchip');
            x_guess=x_guess(:);
            u_guess=interp1(guess.time, guess.inputs,T,'pchip');
            u_guess=u_guess(:);
        else
            x_guess=interp1([guess.time;guess.tf], guess.states,[T;guess.tf],'pchip');
            x_guess=x_guess(:);
            u_guess=interp1(guess.time, guess.inputs,T,'pchip');
            u_guess=u_guess(:);
        end
        data.multipliers=[];
    elseif (strcmp(options.start,'Hot'))
        T=(guess.tf-data.t0)/2*data.tau_inc+(guess.tf+data.t0)/2;
        if size(guess.time,1)==size(guess.states,1)
            x_guess=interp1(guess.time, guess.states,[T;guess.tf],'pchip');
            x_guess=x_guess(:);
            u_guess=interp1(guess.time, guess.inputs,T,'pchip');
            u_guess=u_guess(:);
            data.multipliers.lambda=interp1(guess.time, guess.multipliers.lambda,T,'pchip');
            data.multipliers.lambda=data.multipliers.lambda(:);
            if ng
                lambda_g=interp1(guess.timeFull, guess.multipliers.lambda_g,T,'pchip');
                data.multipliers.lambda=[data.multipliers.lambda;lambda_g(:)];
            end
            if nb
                data.multipliers.lambda=[data.multipliers.lambda;guess.multipliers.lambda_b(1:nb)];
            end
            if nrc
                data.multipliers.lambda=[data.multipliers.lambda;zeros(nrc,1)];
            end
        else
            x_guess=interp1([guess.time;guess.tf], guess.states,[T;guess.tf],'pchip');
            x_guess=x_guess(:);
            u_guess=interp1(guess.time, guess.inputs,T,'pchip');
            u_guess=u_guess(:);
            data.multipliers.lambda=interp1(guess.time, guess.multipliers.lambda,T,'pchip');
            data.multipliers.lambda=data.multipliers.lambda(:);
            if ng
                lambda_g=interp1(guess.time, guess.multipliers.lambda_g,T,'pchip');
                data.multipliers.lambda=[data.multipliers.lambda;lambda_g(:)];
            end
            if nb
                data.multipliers.lambda=[data.multipliers.lambda;guess.multipliers.lambda_b(1:nb)];
            end
            if nrc
                data.multipliers.lambda=[data.multipliers.lambda;zeros(nrc,1)];
            end
        end
    end
else
    if (strcmp(options.start,'Cold'))
        if ~isempty(guess.time)
            Tx=linspace(0,guess.time(end), M);
            Tu=linspace(0,guess.time(end-1), N);
            x_guess=interp1(guess.time, guess.states,Tx,'linear');
            x_guess=reshape(x_guess',M*n,1);
            u_guess=interp1(guess.time, guess.inputs,Tu,'linear');
            u_guess=reshape(u_guess',N*m,1);
        else
            for i=1:n
              x_guess(:,i)=linspace(guess.states(1,i),guess.states(2,i),M);
            end
            for i=1:m
             if N>1  
               u_guess(:,i)=linspace(guess.inputs(1,i),guess.inputs(2,i),N);
              else
               u_guess(:,i)=guess.inputs(1,i);
              end
            end
        end
        data.multipliers=[];
    elseif (strcmp(options.start,'Warm'))
        T=(guess.tf-data.t0)*[0;cumsum(data.tau)]*data.Nm/ns+data.k0;
        if size(guess.time,1)==size(guess.states,1)
            x_guess=interp1(guess.time, guess.states,T,'pchip');
            x_guess=reshape(x_guess',M*n,1);
            u_guess=interp1(guess.time, guess.inputs,T,'pchip');
            u_guess=reshape(u_guess',M*m,1);
        else
            x_guess=interp1([guess.time;guess.tf], guess.states,T,'pchip');
            x_guess=reshape(x_guess',M*n,1);
            u_guess=interp1(guess.time, guess.inputs,T,'pchip');
            u_guess=reshape(u_guess',M*m,1);
        end
        data.multipliers=[];
    elseif (strcmp(options.start,'Hot'))
        T=(guess.tf-data.t0)*[0;cumsum(data.tau)]*data.Nm/ns+data.k0;
        if size(guess.time,1)==size(guess.states,1)
            x_guess=interp1(guess.time, guess.states,T,'pchip');
            x_guess=reshape(x_guess',M*n,1);
            u_guess=interp1(guess.time, guess.inputs,T,'pchip');
            u_guess=reshape(u_guess',M*m,1);
            data.multipliers.lambda=interp1(guess.time, guess.multipliers.lambda,T,'pchip');
            data.multipliers.lambda=reshape(data.multipliers.lambda',M*n,1);
            if ng
                lambda_g=interp1(guess.timeFull, guess.multipliers.lambda_g,T,'linear');
                data.multipliers.lambda=[data.multipliers.lambda;reshape(lambda_g',M*ng,1)];
            end
            if nb
                data.multipliers.lambda=[data.multipliers.lambda;guess.multipliers.lambda_b(1:nb)];
            end
            if nrc
                data.multipliers.lambda=[data.multipliers.lambda;zeros(nrc,1)];
            end
        else
            x_guess=interp1([guess.time;guess.tf], guess.states,T,'pchip');
            x_guess=reshape(x_guess',length(T)*n,1);
            u_guess=interp1(guess.time, guess.inputs,T,'pchip');
            u_guess=reshape(u_guess',length(T)*m,1);
            data.multipliers.lambda=interp1(guess.time, guess.multipliers.lambda,T,'pchip');
            data.multipliers.lambda=reshape(data.multipliers.lambda',M*n,1);
            if ng
                lambda_g=interp1(guess.timeFull, guess.multipliers.lambda_g,T,'linear');
                data.multipliers.lambda=[data.multipliers.lambda;reshape(lambda_g',M*ng,1)];
            end
            if nb
                data.multipliers.lambda=[data.multipliers.lambda;guess.multipliers.lambda_b(1:nb)];
            end
            if nrc
                data.multipliers.lambda=[data.multipliers.lambda;zeros(nrc,1)];
            end
        end
    end
end


if strcmp(options.transcription,'multiple_shooting')
    
%     u_guess=[u_guess;zeros(1,m)];
%     infoNLP.z0=reshape([x_guess u_guess]',M*(n+m),1);
%     infoNLP.z0(end-m+1:end)=[];
%     if nt
%         infoNLP.z0=[guess.tf; infoNLP.z0];
%     end
    if size(u_guess,1)/m==size(x_guess,1)/n
        u_guess=u_guess(1:end-m,:)';u_guess=u_guess(:);
    else
        u_guess=u_guess';u_guess=u_guess(:);
    end
    x_guess=x_guess';x_guess=x_guess(:);
    infoNLP.z0=data.map.xV*x_guess+data.map.uV*u_guess;
else
    u_guess=u_guess';u_guess=u_guess(:);
    x_guess=x_guess';x_guess=x_guess(:);
    infoNLP.z0=data.map.xV*x_guess+data.map.uV*u_guess;
end

if strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')
        if np;infoNLP.z0(n*(M+1)+M*m+1:n*(M+1)+M*m+np)=guess.parameters(:);end
        if nt==1
            infoNLP.z0(nz)=guess.tf;
        elseif nt>=2
            infoNLP.z0(nz-nt+1:nz)=linspace(problem.time.t0,guess.tf,nt);
        end
else
    if np;infoNLP.z0(nt+1:nt+np)=guess.parameters(:);end
    if nt;infoNLP.z0(1)=guess.tf;end
end

%  Formatting options of the ipopt solver

  data.funcs.objective         = @costFunction;
  data.funcs.gradient          = @costGradient;
  data.funcs.constraints       = @constraintFunction;
  data.funcs.jacobian          = @constraintJacobian;
  %data.funcs.iterfunc          = @callback; 
  data.funcs.jacobianstructure = @jacobianstructure;

  
  %Alternative ordering of LGR method 
  if  ((strcmp(options.transcription,'globalLGR')) || (strcmp(options.transcription,'hpLGR'))) && options.reorderLGR
      infoNLP.zl=infoNLP.zl(data.reorder.z_idx);
      infoNLP.zu=infoNLP.zu(data.reorder.z_idx);
      infoNLP.z0=infoNLP.z0(data.reorder.z_idx);
      infoNLP.cl=infoNLP.cl(data.reorder.vert_idx);                    % Lower bounds on constraints.
      infoNLP.cu=infoNLP.cu(data.reorder.vert_idx);                      % Upper bounds on constraints
      data.jacStruct=data.jacStruct(data.reorder.vert_idx,data.reorder.z_idx);
  end
%------------- END OF CODE --------------

