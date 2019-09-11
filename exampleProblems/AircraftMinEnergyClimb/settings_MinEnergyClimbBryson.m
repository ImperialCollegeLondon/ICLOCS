function options = settings_MinEnergyClimbBryson(varargin)

%SETTINGS - General and solver-specific settings are selected here
% Unless specified otherwise the options are set using 0 => no and 1 => yes
%
% Syntax:  options = settings(varargin)
%          When giving one input with varargin, e.g. with settings(20), will use h-method of your choice with N=20 nodes
%          When giving two inputs with varargin, hp-LGR method will be used with two possibilities
%           - Scalar numbers can be used for equally spaced intervals with the same polynoimial orders. For example, settings_hp(5,4) means using 5 LGR intervals each of polynomial degree of 4. 
%           - Alternatively, can supply two arrays in the argument with customized meshing. For example, settings_hp([-1 0.3 0.4 1],[4 5 3]) will have 3 segments on the normalized interval [-1 0.3], [0.3 0.4] and [0.4 1], with polynomial order of 4, 5, and 3 respectively.
%      
% Output:
%    options - Structure containing the settings
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

%------------- BEGIN CODE --------------


%% Transcription Method

% Select a transcription method
%---------------------------------------
% - Direct collocation method ('direct_collocation')
% - Integral residual minimization method ('integral_res_min')
options.transcription='direct_collocation';

% Solution method for integrated residual minimization
%---------------------------------------
% - The alternating method ('alternating')          solving integrated residual minimization and objective minimization in an alternating scheme
% - The penalty barrier method ('weightedCost')     using a weighted objective
options.min_res_mode='alternating';

% For the alternating method, specify priorities for the solution property
%---------------------------------------
% - Lower integrated residual error ('low_res_error')
% - Lower objective value ('low_cost_value')
options.min_res_priority='low_res_error';

% For penalty barrier method, provide a sequence of weights for regularization
%---------------------------------------
options.resCostWeight=[];

% Error criteria (in addition to constraint violation error)
%---------------------------------------
% - local absolute error ('local_abs')      recommended for direct collocation
% - integral residual error ('int_res')     recommended for integral residual minimization
% - both ('both')                           most strict, convergence may be affected (not recommended unless the integrated residual minimization method is selected for result representation )
options.errortype='local_abs';


%% Discretization Method

% Select a discretization method
%---------------------------------------
% Discrete-time model       ('discrete')
% Euler method              ('euler')
% Trapezoidal method        ('trapezoidal')
% Hermite-Simpson method    ('hermite')
% Global LGR method         ('globalLGR')
% Local LGR method          ('hpLGR')
% Automatic chosen direct collocation ('AutoDirect')
options.discretization='hermite';

% Result Representation:
%---------------------------------------
% Direct interpolation in correspondence with the transcription method        ('default')
% Manually selected       ('manual')
% Representation by integrated residual minimization       ('res_min')
% Final representation by integrated residual minimization, with intermediate mesh refinement itrations represented with default method     ('res_min_final_default')      Note: only when error criteria is set to 'both' and the local error tolerance has been fulfilled
% Final representation by integrated residual minimization, with intermediate mesh refinement itrations represented with the manually selected method    ('res_min_final_manual')      Note: only when error criteria is set to 'both' and the local error tolerance has been fulfilled
options.resultRep='default';

% Specify the representation method when 'manual' method is selected
%---------------------------------------
if strcmp(options.resultRep,'manual') || strcmp(options.resultRep,'res_min_final_manual')         
        % State representation
        %---------------------------------------
        %   - Piecewise linear          ('linear'), available for Euler transcription method  
        %   - Piecewise quadratic       ('quadratic'), available for Euler and Trapezoidal transcription methods  
        %   - Piecewise cubic           ('cubic'), available for Hermite-Simpson transcription method  
        %   - Barycentric Lagrange Interpolation ('Barycentric'), available for LGR transcription method  
        %   - Legendre polynomial fitting  ('Legendre'), available for LGR transcription method  
        %   - Piecewise Cubic Hermite Interpolating Polynomial with Matlab pchip function        ('pchip'), available for all transcription methods
        options.stateRep='pchip';
        % Input representation
        %---------------------------------------
        %   - Piecewise constant        ('constant'), available for all transcription methods
        %   - Piecewise linear          ('linear'), available for all transcription methods
        %   - Piecewise quadratic       ('quadratic'), available for Trapezoidal transcription methods  
        %   - Barycentric Lagrange Interpolation ('Barycentric'), available for LGR transcription method  
        %   - Legendre polynomial fitting  ('Legendre'), available for LGR transcription method  
        %   - Piecewise Cubic Hermite Interpolating Polynomial with Matlab pchip function        ('pchip'), available for all transcription methods
        options.inputRep='linear';
end

% Further settings for solution representation by integrated residual minimization
%---------------------------------------
if strcmp(options.resultRep,'res_min') || strcmp(options.resultRep,'res_min_final_default') || strcmp(options.resultRep,'res_min_final_manual')
        % Option on whether to match the end-point results from collocation, recommend to select No (0) 
        options.resminRep.collmatch=0; 
        % Tolerance for the cost value, in precentage (e.g. 0.01 represent maximum increase of 1% of the original cost value)
        options.resminRep.costTol=0; %
end

%% Derivative generation

% Derivative computation method
%---------------------------------------
% Analytic differentiation: analytic gradients   ('analytic')
    % Whenever the analytic differentiation is enabled it is necessary to specify the available analytic forms for the cost function, the dynamic equations and the constraints in the appropriate files .m
% Numerical differentiation: finite differences  ('numeric')
% Algorithmic differentiation with Adigator      ('adigator')
    % Make sure you provide the path to the Adigator directory of startupadigator.m
options.derivatives='numeric';
options.adigatorPath='../../adigator';

% Perturbation sizes for numerical differentiation
%---------------------------------------
%  It is possible to select default values for the perturbations by setting  options.perturbation.H and 
%  options.perturbation.J to the empty matrix.
%  The default values for the gradient approximation is (eps/2)^(1/3)
%  while for the  second derivative is (8*eps)^(1/3). 
options.perturbation.H=[];  % Perturbation size for the second derivatives
options.perturbation.J=[];  % Perturbation size for the first derivatives

%% NLP solver

% Select a NLP solver
%---------------------------------------
% IPOPT: recommended but needs ipopt.mex        ('ipopt')
% fmincon                                       ('fmincon')
% WORHP                                         ('worhp')
options.NLPsolver='ipopt';

% IPOPT settings (if required)
%---------------------------------------
options.ipopt.tol=1e-9;                        % Desired convergence tolerance (relative). The default value is  1e-8. 
options.ipopt.print_level=5;                   % Print level. The valid range for this integer option is [0,12] and its default value is 5.
options.ipopt.max_iter=5000;                   % Maximum number of iterations. The default value is 3000.
 
options.ipopt.mu_strategy ='adaptive';         % Determines which barrier parameter update strategy is to be used. 
                                               % The default value for this string option is "monotone".
                                               % Possible values:
                                               %   'monotone': use the monotone (Fiacco-McCormick) strategy
                                               %   'adaptive': use the adaptive update strategy

options.ipopt.hessian_approximation='exact';   %  Indicates what information for the Hessian of the Lagrangian function is                                                    
                                               %  used by the algorithm. The default value is 'exact'.
                                               %  Possible values:
                                               %   'exact': Use second derivatives provided by ICLOCS.
                                               %   'limited-memory': Perform a limited-memory quasi-Newton approximation
					                           %		             implemented inside IPOPT

options.ipopt.limited_memory_max_history=6;   % Maximum size of the history for the limited quasi-Newton Hessian approximation. The valid range for this integer option is [0, +inf) 
                                               % and its default value is 6. 
options.ipopt.limited_memory_max_skipping=1;  % Threshold for successive iterations where update is skipped for the quasi-Newton approximation.
                                               % The valid range for this integer option is [1,+inf) and its default value is 2. 

% fmincon settings (NOT RECOMMENDED!)
%---------------------------------------
% See website for detailed info

% WORHP settings
%---------------------------------------
% need to be configured with the xml file


%% Meshing Strategy

% Type of meshing
%---------------------------------------
% - fixed mesh ('fixed')
% - with local refinement of mesh ('mesh_refinement')
% - flexible mesh with adaptively spaced segments, ONLY AVAILABLE FOR hpLGR Discretization with Direct Collocation Method ('hp_flexible')
options.meshstrategy='mesh_refinement';

% Mesh Refinement Method
%---------------------------------------
% Increase Polynomial Order        ('IO')
% Add intervals                    ('AI')
% Automatic refinement             ('Auto')
% For current version leave the setting to Automatic
options.MeshRefinement='Auto';

% Mesh Refinement Preferences
%---------------------------------------
% Prioritize MR time               ('aggressive')   A relative aggressive scheme that aim to reduce the number of MR iterations and the total MR time
% Prioritize MR efficiency         ('efficient')   A relative aggressive scheme that aim to reduce the size of the problem at the end of MR iterations, making it potentially more efficient for online re-computations
options.MRstrategy='aggressive';

% Maximum number of mesh refinement iterations
%---------------------------------------
options.maxMRiter=50;

% Discountious Input
%---------------------------------------
% NOTE: Not yet available in this version
options.disContInputs=0;

% Minimum and maximum time interval
%---------------------------------------
% Define the minimum and maximum time interval (in the same unit as in t) for the hp-flexible method and for mesh refinement schemes
options.mintimeinterval=1e-09; 
options.maxtimeinterval=inf; 

% Distribution of integration steps. 
%---------------------------------------
% Set tau=0 for equispaced steps.
% Otherwise: tau is a vector of length M-1 with 0<tau(i)<1 and sum(tau)=1.
% For discrete time system  set  tau=0.
options.tau=0;


%% Other Settings

% Cold/Warm/Hot Start 
%---------------------------------------
options.start='Cold';

% Automatic scaling (recommended)
%---------------------------------------
options.scaling=1;

% Reorder of LGR Method
%---------------------------------------
% NOTE: For current version, leave it to off (0)
options.reorderLGR=0;

% Early termination of residual minimization if tolerance is met
%---------------------------------------
% For use with integral residual minimization method with the alternating scheme ONLY
options.resminEarlyStop=0;

% External Constraint Handling
%---------------------------------------
% For problems with large sets of inactive constraints, external constraint handling may help to reduce the compuatation cost
% Enable or disable the external constraint handling function (1 or 0)
options.ECH.enabled=0; 
% A conservative buffer zone (default to 0.1, i.e. 10%)
options.ECH.buffer_pct=0.1; 

% Regularization Strategy
%---------------------------------------
% A regularization parameter could be setup in the problem file with problem.data.penalty.values defining the regularization sequence of weights, and problem.data.penalty.i the starting index. The parameter could be called with vdat.penalty.values(vdat.penalty.i) in problem formulation.
% Here, multiple strategies are possible interacting with mesh refinement
% - Off ('off')                                 
% - Regularization priority ('reg_priority')    First change regularization paramters iteratively with warm starting resolve, then perform mesh refinement iterations
% - Mesh Refinement priority ('MR_priority')    First perform mesh refinement iterations, then change regularization paramters iteratively with warm starting resolve
% - Simultaneous ('simultaneous')               Perform mesh refinement and regulzarization paramter iterations simultaneously
options.regstrategy='off';

% Auto selection of h/hp method based on input formulation of the settings function call
%---------------------------------------
% LEAVE THIS PART UNCHANGED AND USE FUNCTION SYNTAX (AS DESCRIBED ON THE TOP) TO DEFINE THE ITEGRATION NODES
if nargin==2
    if strcmp(varargin{2},'h')
        options.nodes=varargin{1}; 
        options.discretization='hermite';
    else
        if length(varargin{1})==1
            options.nsegment=varargin{1}; 
            options.pdegree=varargin{2}; 
        else
            options.tau_segment=varargin{1};
            options.npsegment=varargin{2};
        end
        options.discretization='hpLGR';
    end
else
    options.nodes=varargin{1}; 
end


%% Output settings

% Display computation time
%---------------------------------------
options.print.time=1;

% Display relative local discretization error (recommended for direct transcription)
%---------------------------------------
options.print.relative_local_error=1;

% Display cost (objective) values
%---------------------------------------
options.print.cost=1;


% Plot figures £¨when genSolutionPlots function is called)
%---------------------------------------
% 0: Do not plot
% 1: Plot all figures (state and input trajectory, multipliers/costate values and errors)
% 2: Plot only the state and input trajectory
% 3: Plot only the multipliers/costate values
% 4: Plot only the error values (absolute local error, relative local error and absolute constraint violation error)
options.plot=1;


