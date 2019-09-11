function config_Check(varargin)
%config_Check - check for some incompatible settings
%
% Syntax:   config_Check(options,nps,npd)  (hp method)
%           config_Check(options,N)        (h method)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

if nargin==3
    options=varargin{1};
    nps=varargin{2};
    npd=varargin{3};
    
    if strcmp(options.transcription,'integral_res_min')
        if strcmp(options.derivatives,'adigator')
            if nps<2 || any(npd<2)
                error('Configuration not suitable for integrated residual minimization with hp method. Need atleast 2 mesh segments each with minimum polynomial order of 2.')
            end
        elseif strcmp(options.derivatives,'numeric')
            warning('Numeric differentiation can be slow for integrated residual minimization, please consider to use algorithmic differentiation toolbox of Adigator')
            if nps<3 || any(npd<3)
                error('Configuration not suitable for integrated residual minimization with hp method. Need atleast 2 mesh segments each with minimum polynomial order of 3 (2 segments and order of 2 incase of using Adigator).')
            end
        else
            error('Analytic derivatives not yet supported for integrated residual minimization, please consider to use algorithmic differentiation toolbox of Adigator.')
        end
    elseif strcmp(options.transcription,'direct_collocation')
        if any(npd<2)
            error('Do not supported polynomial order less than 2.')
        end
    end
else
    options=varargin{1};
    N=varargin{2};
    if N<3
        error('Do not supported number of mesh nodes less than 3.')
    end
    if strcmp(options.transcription,'integral_res_min')
        if strcmp(options.derivatives,'numeric')
            warning('Numeric differentiation can be slow for integrated residual minimization, please consider to use algorithmic differentiation toolbox of Adigator')
        elseif strcmp(options.derivatives,'analytic')
            error('Analytic derivatives not yet supported for integrated residual minimization, please consider to use algorithmic differentiation toolbox of Adigator.')
        end
    elseif strcmp(options.transcription,'direct_collocation')

    end
end

if strcmp(options.derivatives,'adigator') && ~strcmp(options.NLPsolver,'ipopt') 
    error('Please use the NLP solver of IPOPT for derivative information supplied with Adigator toolbox')
end

end

