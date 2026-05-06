function tests = testInstaller
%TESTINSTALLER Smoke tests for one-command ICLOCS setup.
tests = functiontests(localfunctions);
end

function testInstallIclocsConfiguresPathAndSmokeCheck(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
cleanupRoot = onCleanup(@() removePath(rootDir));

report = installIclocs('savePath',false,'runSmokeTest',true,'quiet',true);
cleanupInstall = onCleanup(@() removeInstallPaths(rootDir,report));

verifyEqual(testCase,report.repoRoot,rootDir);
verifyEqual(testCase,exist('solveMyProblem','file'),2);
verifyTrue(testCase,report.ready);
verifyTrue(testCase,report.smokeCheck.passed);

delete(cleanupInstall);
delete(cleanupRoot);
end

function testIclocsDoctorReportsDependencyState(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
toolsDir = fullfile(rootDir,'tools');
addpath(toolsDir);
cleanupTools = onCleanup(@() removePath(toolsDir));

setupIclocsPath(rootDir);
cleanupSource = onCleanup(@() removePath(genpath(fullfile(rootDir,'src'))));

report = iclocsDoctor(rootDir,'quiet',true);

verifyEqual(testCase,report.repoRoot,rootDir);
verifyTrue(testCase,isfield(report,'dependencies'));
verifyTrue(testCase,numel(report.dependencies) >= 5);
verifyTrue(testCase,isfield(report,'solverAvailable'));
verifyTrue(testCase,report.ready);

delete(cleanupSource);
delete(cleanupTools);
end

function removeInstallPaths(rootDir,report)
removePath(genpath(fullfile(rootDir,'src')));
removePath(fullfile(rootDir,'tools'));
if isfield(report,'examplePath') && ~isempty(report.examplePath)
    removePath(report.examplePath);
end
removePath(rootDir);
end

function removePath(pathToRemove)
try
    rmpath(pathToRemove);
catch
end
end
