function report = installIclocs(varargin)
%INSTALLICLOCS Configure ICLOCS and check local MATLAB dependencies.
%
% Syntax:
%   installIclocs
%   installIclocs('savePath',true)
%   installIclocs('runSmokeTest',false)
%   report = installIclocs(...)
%
% Name-value inputs:
%   savePath     - Save the updated MATLAB path with SAVEPATH. The default
%                  is false so first-time runs do not modify user startup
%                  state unless requested.
%   runSmokeTest - Build a small example problem without solving it. The
%                  default is true.
%   example      - Example problem used for the smoke check. The default is
%                  'DoubleIntegratorTracking'.
%   quiet        - Suppress progress output. The default is false.
%
% Outputs:
%   report - Structure containing path, dependency, and smoke-check status.

parser = inputParser;
parser.FunctionName = 'installIclocs';
addParameter(parser,'savePath',false,@isFlag);
addParameter(parser,'runSmokeTest',true,@isFlag);
addParameter(parser,'example','DoubleIntegratorTracking',@isTextScalar);
addParameter(parser,'quiet',false,@isFlag);
parse(parser,varargin{:});

options = parser.Results;
options.savePath = logical(options.savePath);
options.runSmokeTest = logical(options.runSmokeTest);
options.example = char(options.example);
options.quiet = logical(options.quiet);

repoRoot = fileparts(mfilename('fullpath'));
toolsDir = fullfile(repoRoot,'tools');
addpath(toolsDir);

setupIclocsPath(repoRoot);
report = iclocsDoctor(repoRoot,'quiet',true);
report.example = options.example;

exampleDir = fullfile(repoRoot,'exampleProblems',options.example);
if exist(exampleDir,'dir') == 7
    addpath(exampleDir);
    report.examplePath = exampleDir;
else
    report.examplePath = '';
    report.warnings{end+1} = ['Example folder not found: ' exampleDir];
end

if options.runSmokeTest
    report.smokeCheck = runSmokeCheck(options.example, exampleDir);
else
    report.smokeCheck = struct('example',options.example,'passed',true, ...
        'message','Skipped by request.');
end

report.savePath = struct('requested',options.savePath,'status',[],'message','');
if options.savePath
    [status,message] = savepath;
    report.savePath.status = status;
    report.savePath.message = message;
    if status ~= 0
        report.warnings{end+1} = ['Could not save MATLAB path: ' message];
    end
end

if ~options.quiet
    printInstallSummary(report);
end

if nargout == 0
    clear report
end
end

function smokeCheck = runSmokeCheck(exampleName, exampleDir)
smokeCheck = struct('example',exampleName,'passed',false,'message','', ...
    'problemFields',{{}},'guessFields',{{}});

if exist(exampleDir,'dir') ~= 7
    smokeCheck.message = ['Example folder not found: ' exampleDir];
    return
end

try
    [problem,guess] = feval(exampleName);
    requiredProblemFields = {'settings','time','states','inputs','constraints'};
    requiredGuessFields = {'time','states','inputs'};
    smokeCheck.problemFields = requiredProblemFields;
    smokeCheck.guessFields = requiredGuessFields;
    smokeCheck.passed = hasFields(problem,requiredProblemFields) && ...
        hasFields(guess,requiredGuessFields);
    if smokeCheck.passed
        smokeCheck.message = 'Problem and guess structures were created.';
    else
        smokeCheck.message = 'Problem or guess structure is missing expected fields.';
    end
catch err
    smokeCheck.message = err.message;
end
end

function tf = hasFields(value,fieldNames)
tf = isstruct(value);
idx = 1;
while tf && idx <= numel(fieldNames)
    tf = isfield(value,fieldNames{idx});
    idx = idx + 1;
end
end

function printInstallSummary(report)
fprintf('\nICLOCS setup\n');
fprintf('  Repository: %s\n',report.repoRoot);
fprintf('  Source path: %s\n',statusWord(report.srcOnPath));
fprintf('  Tools path: %s\n',statusWord(report.toolsOnPath));
fprintf('  MATLAB version: %s\n',statusWord(report.matlabOk));
fprintf('  Solver available: %s\n',statusWord(report.solverAvailable));
fprintf('  Smoke check: %s\n',statusWord(report.smokeCheck.passed));

if ~isempty(report.warnings)
    fprintf('\nWarnings\n');
    for idx = 1:numel(report.warnings)
        fprintf('  - %s\n',report.warnings{idx});
    end
end

fprintf('\nNext step\n');
if report.solverAvailable
    fprintf('  cd(''%s'')\n',report.examplePath);
    fprintf('  main_%s\n',report.example);
else
    fprintf('  Install IPOPT, Optimization Toolbox, or WORHP before solving examples.\n');
    fprintf('  See docs/INSTALLATION.md for solver setup notes.\n');
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

function tf = isFlag(value)
tf = (islogical(value) || isnumeric(value)) && isscalar(value);
end

function tf = isTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end
