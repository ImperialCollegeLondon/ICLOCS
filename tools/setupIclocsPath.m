function repoRoot = setupIclocsPath(repoRoot)
%SETUPICLOCSPATH Add ICLOCS source folders to the MATLAB path.
%
% Syntax:
%   setupIclocsPath
%   setupIclocsPath(repoRoot)
%   repoRoot = setupIclocsPath(...)
%
% Inputs:
%   repoRoot - Optional path to the ICLOCS repository root. When omitted,
%              the path is inferred from this file location.
%
% Outputs:
%   repoRoot - Repository root that was added.

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end

addpath(genpath(fullfile(repoRoot,'src')));

if nargout == 0
    clear repoRoot
end
end
