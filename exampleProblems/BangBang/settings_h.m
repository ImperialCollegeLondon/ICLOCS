function options = settings_h(N)

%SETTINGS - General and solver-specific settings are selected here
% Unless specified otherwise the options are set using 0 => no and 1 => yes
%
% Syntax:  options = settings
%
% Output:
%    options - Structure containing the settings
%
% Other m-files required: none
% Subfunctions: none
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

%------------- BEGIN CODE --------------


% Transcription Method:
%---------------------------------------
% Discrete-time model       ('discrete')
% Multiple shooting method  ('multiple_shooting') WARNING: The
%                                                 'quasi-newton' option for the hessian 
%						                           computation has to be selected
%                                                 (options.ipopt.hessian_approximation='limited-memory').
% Euler method              ('euler')
% Trapezoidal method        ('trapezoidal')
% Hermite-Simpson method    ('hermite')
% Global LGR method         ('globalLGR')
% Local LGR method          ('hpLGR')
% Automatic chosen direct collocation ('AutoDirect')
options.transcription='hermite';

% Result Representation:
%---------------------------------------
% As recommended        ('default')
% Direct reconstrunction in correspondence with the transcription method       ('direct')
% Manually select       ('manual')
options.resultRep='default';

% Maunal selection of result representation method:
%---------------------------------------
% State representation
%   - Piecewise linear          ('linear'), available for Euler transcription method  
%   - Piecewise quadratic       ('quadratic'), available for Euler and Trapezoidal transcription methods  
%   - Piecewise cubic           ('cubic'), available for Hermite-Simpson transcription method  
%   - Legendre polynomials      ('Legendre'), available for LGR transcription method  
%   - Piecewise Cubic Hermite Interpolating Polynomial with Matlab pchip function        ('pchip'), available for all transcription methods
options.stateRep='cubic';
% Input representation
%   - Piecewise constant        ('constant'), available for all transcription methods
%   - Piecewise linear          ('linear'), available for all transcription methods
%   - Piecewise quadratic       ('quadratic'), available for Trapezoidal transcription methods  
%   - Legendre polynomials      ('Legendre'), available for LGR transcription method  
%   - Piecewise Cubic Hermite Interpolating Polynomial with Matlab pchip function        ('pchip'), available for all transcription methods
options.inputRep='quadratic';

% Derivative generation :
%---------------------------------------
% Whenever the analytic differentiation is enabled it is necessary to
% specify the available analytic forms for the cost function, the dynamic equations 
% and  the constraints in the appropriate files .m

% Numerical differentiation: finite differences  ('numeric')
% Analytic differentiation: analytic gradients   ('analytic')
% Algorithmic differentiation with adigator  ('adigator'). Make sure you
% run startupadigator.m first in the adigator directory every time the MATLAB is restarted.
options.derivatives='numeric';


% Numeric generation of the Hessian:
%----------------------------------------------------------------

% Whenever the numeric differentiation is enabled it is necessary to
% specify which kind of finite difference approximation to use  between 
% the following ones:
% 
% Central difference ('central')  
options.hessianFD='central';


%  The perturbation size for numerical second derivatives 
%  can be set in options.perturbation.H. The perturbation size for numerical first derivatives 
%  can be set in  options.perturbation.J. 
%  It is possible to select default values for the perturbations by setting  options.perturbation.H and 
%  options.perturbation.J to the empty matrix.
%  The default values for the gradient approximation is (eps/2)^(1/3)
%  while for the  second derivative is (8*eps)^(1/3). 

options.perturbation.H=[];  % Perturbation size for the second derivatives
options.perturbation.J=[];  % Perturbation size for the first derivatives

% NLP solver
%---------------------------------------
% IPOPT: recommended but needs ipopt.mex        ('ipopt')
% fmincon                                       ('fmincon')
% WORHP                                         ('worhp')
options.NLPsolver='ipopt';

% IPOPT settings (if required)
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
% See website for detailed info

% WORHP settings needed to be configured with the xml file

% Cold/Warm/Hot Start (recommended)
%---------------------------------------
options.start='Cold';

% Automatic scaling (recommended)
%---------------------------------------
options.scaling=1;


% Output settings
%---------------------------------------

% Display computation time
options.print.time=1;

% Display relative local discretization error (recommended for direct transcription)
options.print.relative_local_error=1;

% Display cost
options.print.cost=1;


% Plot states
options.plot.states=1;

% Plot inputs
options.plot.inputs=1;

% Plot Lagrange multipliers
options.plot.multipliers=1;


% Direct transcription settings
%---------------------------------------


% Number of   integration nodes in the interval t=[0,tf]; nodes=steps+1.
% The quantity steps/N (N number of control actions) must be a positive
% integer. For LGR: Number of LGR points on interval t=[-1,tau_n], tau_n<1
options.nodes=N; 

% Adaptively spaced segments
options.adaptseg=0; 

% Minimum time interval
options.mintimeinterval=1e-09; 


% Distribution of integration steps. Set tau=0 for equispaced steps.
% Otherwise: tau is a vector of length M-1 with 0<tau(i)<1 and sum(tau)=1.
% For discrete time system  set  tau=0.

options.tau=0;


% Multiple shooting settings
%---------------------------------------
% N/S=normal/stiff. H/M/L=high/medium/low accuracy
%
% 'cvodes'   N/S    H    High accuracy, difficult problems(slow)
%
%
%  Note: 'cvodes' requires the sundialsTB. 
%        

% 
% options.ODEsolver='cvodes';
% 
% 
% 
% 
% % CVODES settings (if required)
% 
% Method='Adams';  % Method: Adams, BDF
% Solver='Newton'; % Solver: Newton, Functional (requires dfdx)
% 
% 
% % Options for forward integratio nwhen cvodes is enabled
% options.cvodes = CVodeSetOptions('RelTol',1.e-4,...
%                           'AbsTol',1.e-6,...
%                           'LinearSolver','Dense',...
%                           'MaxNumSteps',10000,...
%                           'LMM',Method,...
%                           'NonlinearSolver',Solver); 
%                  
% % Forward sensitivity options when cvodes is enabled 
% 
% options.cvodesf=CVodeSensSetOptions('ErrControl',true,'method','Staggered'); % FSA initialization
               

















