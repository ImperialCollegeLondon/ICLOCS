function [ ] = genAdigator4ICLOCS( options, probinfo, n, m, np )
%genAdigator4ICLOCS - Initialize Adigator for use with ICLOCS
%
% Syntax:  [ ] = genAdigator4ICLOCS( options, probinfo, n, m, nt, np )
%
% Inputs:
%    required - Flag that determines what to compute for current z
%    z - Current NLP variable
%    data - Structure of data required to compute
%
% Outputs:
%    solution - Data structure containing the solution
%
% Other m-files required: jacobianFD_LGR, hessianCD_LGR (optional)
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

% ------------------------------ ADiGator ------------------------------- %
gX = adigatorCreateDerivInput([Inf,n],struct('vodname','Y','vodsize',[Inf,(m+n)],...
  'nzlocs',[(1:n).' (1:n).']));
gU = adigatorCreateDerivInput([Inf,m],struct('vodname','Y','vodsize',[Inf,(m+n)],...
  'nzlocs',[(1:m).' (n+1:n+m).']));
gT = adigatorCreateDerivInput([Inf,1],struct('vodname','T','vodsize',[Inf,1],...
  'nzlocs',[(1:1).' (1:1).']));
if np
   gP = adigatorCreateDerivInput([Inf,n],struct('vodname','P','vodsize',[Inf,(m+n)],...
  'nzlocs',[(1:np).' (1:np).']));
else
   gP=[];
end

% Create Vectorized Dynamics 1st Derivative File
adigator(probinfo.data.plantmodel,{gX,gU,gP,gT,probinfo.data},'dynamics_Y',adigatorOptions('overwrite',1));


if ~strcmp(options.ipopt.hessian_approximation,'limited-memory')
    % Create Vectorized Dynamics 2nd Derivative File
    gXs.f = gX; gUs.f = gU; gTs.f=gT; gPs.f=gP;
    % Need to make aux inputs for derivative fields since we dont know their
    % vectorized dimension
    gXs.dY = adigatorCreateAuxInput([Inf n]);
    gUs.dY = adigatorCreateAuxInput([Inf m]);
    gTs.dT = adigatorCreateAuxInput([Inf 1]);
    if np
        gPs.dY = adigatorCreateAuxInput([Inf np]);
    else
        gPs=[];
    end
    adigator('dynamics_Y',{gXs,gUs,gPs,gTs,probinfo.data},'dynamics_YY',...
      adigatorOptions('overwrite',1,'comments',0));
end

end

