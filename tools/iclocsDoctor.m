function report = iclocsDoctor(repoRoot,varargin)
%ICLOCSDOCTOR Report ICLOCS path and dependency status.
%
% Syntax:
%   iclocsDoctor
%   iclocsDoctor(repoRoot)
%   report = iclocsDoctor(...,'quiet',true)
%
% Inputs:
%   repoRoot - Optional path to the ICLOCS repository root. When omitted,
%              the path is inferred from this file location.
%
% Name-value inputs:
%   quiet - Suppress printed output. The default is false.
%
% Outputs:
%   report - Structure with MATLAB version, path, optional dependency, and
%            solver availability checks.

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
elseif isParameterName(repoRoot)
    varargin = [{repoRoot} varargin];
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end

parser = inputParser;
parser.FunctionName = 'iclocsDoctor';
addParameter(parser,'quiet',false,@isFlag);
parse(parser,varargin{:});
quiet = logical(parser.Results.quiet);

report = struct;
report.repoRoot = char(repoRoot);
report.matlabVersion = version;
report.matlabRelease = matlabReleaseName;
report.matlabOk = matlabVersionAtLeast('9.2');
report.toolsOnPath = exist('setupIclocsPath','file') == 2;
report.srcOnPath = exist('solveMyProblem','file') == 2;
report.dependencies = dependencyChecks(report.repoRoot);
report.solverAvailable = hasAnySolver(report.dependencies);
report.ready = report.matlabOk && report.toolsOnPath && report.srcOnPath;
report.readyToSolve = report.ready && report.solverAvailable;
report.warnings = {};

if ~report.matlabOk
    report.warnings{end+1} = 'MATLAB R2017a or later is recommended for ICLOCS 2.5.';
end

if ~report.srcOnPath
    report.warnings{end+1} = 'ICLOCS source folders are not on the MATLAB path. Run setupIclocsPath.';
end

if ~report.solverAvailable
    report.warnings{end+1} = 'No supported NLP solver was detected. Install IPOPT, Optimization Toolbox, or WORHP before solving examples.';
end

if ~quiet
    printDoctorReport(report);
end

if nargout == 0
    clear report
end
end

function dependencies = dependencyChecks(repoRoot)
dependencies = emptyDependency;

dependencies(end+1) = makeDependency( ...
    'IPOPT MATLAB interface', ...
    hasFunctionOrMex('ipopt'), ...
    false, ...
    'Preferred NLP solver for most ICLOCS examples.', ...
    'Install a MATLAB-compatible IPOPT interface and add it to the MATLAB path.');

dependencies(end+1) = makeDependency( ...
    'Optimization Toolbox', ...
    hasFunctionOrMex('fmincon') && hasFunctionOrMex('optimoptions'), ...
    false, ...
    'Provides fmincon for small sanity-check solves.', ...
    'Install MathWorks Optimization Toolbox or choose another supported solver.');

dependencies(end+1) = makeDependency( ...
    'ADiGator', ...
    hasFunctionOrMex('adigator') || exist(fullfile(repoRoot,'adigator'),'dir') == 7, ...
    false, ...
    'Optional algorithmic differentiation dependency.', ...
    'Install ADiGator and set options.adigatorPath when using ADiGator derivatives.');

dependencies(end+1) = makeDependency( ...
    'WORHP MATLAB interface', ...
    hasFunctionOrMex('worhp') || hasFunctionOrMex('worhp_interface'), ...
    false, ...
    'Optional supported NLP solver.', ...
    'Install a WORHP release with the MATLAB interface if your workflow uses WORHP.');

dependencies(end+1) = makeDependency( ...
    'Simulink', ...
    hasFunctionOrMex('simulink') || safeLicenseTest('Simulink'), ...
    false, ...
    'Optional dependency for workflows that use Simulink models.', ...
    'Install Simulink only if your model or example requires it.');
end

function dependency = emptyDependency
dependency = struct('name',{},'available',{},'required',{},'detail',{},'help',{});
end

function dependency = makeDependency(name,available,required,detail,helpText)
dependency = struct('name',name,'available',logical(available), ...
    'required',logical(required),'detail',detail,'help',helpText);
end

function tf = hasAnySolver(dependencies)
solverNames = {'IPOPT MATLAB interface','Optimization Toolbox','WORHP MATLAB interface'};
tf = false;
idx = 1;
while ~tf && idx <= numel(dependencies)
    tf = dependencies(idx).available && any(strcmp(dependencies(idx).name,solverNames));
    idx = idx + 1;
end
end

function tf = hasFunctionOrMex(functionName)
location = exist(functionName,'file');
tf = location == 2 || location == 3 || location == 6;
end

function tf = safeLicenseTest(featureName)
try
    tf = license('test',featureName);
catch
    tf = false;
end
tf = logical(tf);
end

function tf = matlabVersionAtLeast(versionNumber)
try
    tf = ~verLessThan('matlab',versionNumber);
catch
    tf = true;
end
end

function releaseName = matlabReleaseName
try
    releaseName = version('-release');
catch
    releaseName = '';
end
end

function printDoctorReport(report)
fprintf('\nICLOCS doctor\n');
fprintf('  Repository: %s\n',report.repoRoot);
fprintf('  MATLAB: %s %s [%s]\n',report.matlabVersion,report.matlabRelease,statusWord(report.matlabOk));
fprintf('  tools on path: %s\n',statusWord(report.toolsOnPath));
fprintf('  src on path: %s\n',statusWord(report.srcOnPath));

fprintf('\nDependencies\n');
for idx = 1:numel(report.dependencies)
    dependency = report.dependencies(idx);
    fprintf('  [%s] %s - %s\n',statusWord(dependency.available), ...
        dependency.name,dependency.detail);
    if ~dependency.available
        fprintf('       %s\n',dependency.help);
    end
end

if ~isempty(report.warnings)
    fprintf('\nWarnings\n');
    for idx = 1:numel(report.warnings)
        fprintf('  - %s\n',report.warnings{idx});
    end
end
fprintf('\n');
end

function word = statusWord(value)
if value
    word = 'ok';
else
    word = 'missing';
end
end

function tf = isParameterName(value)
tf = false;
if ischar(value)
    tf = strcmp(value,'quiet');
elseif isstring(value) && isscalar(value)
    tf = strcmp(char(value),'quiet');
end
end

function tf = isFlag(value)
tf = (islogical(value) || isnumeric(value)) && isscalar(value);
end
