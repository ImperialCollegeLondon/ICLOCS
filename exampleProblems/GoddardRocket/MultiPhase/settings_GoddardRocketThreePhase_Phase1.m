function options = settings_GoddardRocketThreePhase_Phase1(varargin)

%SETTINGS - General and solver-specific settings are selected here (phase specific)
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



% Discretization Method:
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
% Representation by integrated residual minimization       ('res_min')
% Representation by integrated residual minimization, only when error criteria is set to 'both' and the local error tolerance has been fulfilled    ('res_min_final_default')
% Representation by integrated residual minimization, only when error criteria is set to 'both' and the local error tolerance has been fulfilled    ('res_min_final_manual')
% Manually selected       ('manual')
options.resultRep='default';

% Further settings
if strcmp(options.resultRep,'res_min') || strcmp(options.resultRep,'res_min_final_default') || strcmp(options.resultRep,'res_min_final_manual')
        % Matching the end-point results from collocation (1)
        % Do not match results from collocation, i.e. only use as intial values (0)
        options.resminRep.collmatch=0; 
        options.resminRep.costTol=0; %
end
if strcmp(options.resultRep,'manual') || strcmp(options.resultRep,'res_min_final_manual') 
        % Maunal selection of result representation method:
        %---------------------------------------
        % State representation
        %   - Piecewise linear          ('linear'), available for Euler transcription method  
        %   - Piecewise quadratic       ('quadratic'), available for Euler and Trapezoidal transcription methods  
        %   - Piecewise cubic           ('cubic'), available for Hermite-Simpson transcription method  
        %   - Barycentric Lagrange Interpolation ('Barycentric'), available for LGR transcription method  
        %   - Legendre polynomial fitting  ('Legendre'), available for LGR transcription method  
        %   - Piecewise Cubic Hermite Interpolating Polynomial with Matlab pchip function        ('pchip'), available for all transcription methods
        options.stateRep='pchip';
        % Input representation
        %   - Piecewise constant        ('constant'), available for all transcription methods
        %   - Piecewise linear          ('linear'), available for all transcription methods
        %   - Piecewise quadratic       ('quadratic'), available for Trapezoidal transcription methods  
        %   - Barycentric Lagrange Interpolation ('Barycentric'), available for LGR transcription method  
        %   - Legendre polynomial fitting  ('Legendre'), available for LGR transcription method  
        %   - Piecewise Cubic Hermite Interpolating Polynomial with Matlab pchip function        ('pchip'), available for all transcription methods
        options.inputRep='linear';
end


% Direct transcription settings
%---------------------------------------


% Number of integration nodes in the interval t=[0,tf]; nodes=steps+1.
% The quantity steps/N (N number of control actions) must be a positive
% integer. For LGR: Number of LGR points on interval t=[-1,tau_n], tau_n<1
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

% Adaptively spaced segments
options.adaptseg=0; 

% Minimum time interval
options.mintimeinterval=0.1; 


% Distribution of integration steps. Set tau=0 for equispaced steps.
% Otherwise: tau is a vector of length M-1 with 0<tau(i)<1 and sum(tau)=1.
% For discrete time system  set  tau=0.

options.tau=0;

