function tests = testSetupIclocsPath
%TESTSETUPICLOCSPATH Smoke test for the repository path setup helper.
tests = functiontests(localfunctions);
end

function testSetupIclocsPathAddsSourceTree(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
toolsDir = fullfile(rootDir,'tools');
addpath(toolsDir);
cleanupTools = onCleanup(@() rmpath(toolsDir));

repoRoot = setupIclocsPath(rootDir);
cleanupSource = onCleanup(@() rmpath(genpath(fullfile(rootDir,'src'))));

verifyEqual(testCase, repoRoot, rootDir);
verifyEqual(testCase, exist('solveMyProblem','file'), 2);

delete(cleanupSource);
delete(cleanupTools);
end
