function tests = testExamples
%TESTEXAMPLES Smoke tests for example definitions without running solvers.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(rootDir,'src')));
testCase.TestData.RootDir = rootDir;
end

function testDoubleIntegratorTrackingDefinition(testCase)
exampleDir = fullfile(testCase.TestData.RootDir,'exampleProblems','DoubleIntegratorTracking');
cleanup = addExamplePath(exampleDir);

[problem,guess] = DoubleIntegratorTracking;

verifyEqual(testCase, problem.time.t0_min, 0);
verifyEqual(testCase, problem.time.tf_min, 10);
verifyEqual(testCase, problem.time.tf_max, 10);
verifyEqual(testCase, problem.states.x0, [0 5]);
verifyEqual(testCase, problem.inputs.ul, -10);
verifyEqual(testCase, problem.inputs.uu, 10);
verifySize(testCase, guess.states, [5 2]);
verifySize(testCase, guess.inputs, [5 1]);
verifyTrue(testCase, isa(problem.settings,'function_handle'));
verifyTrue(testCase, isa(problem.data.InternalDynamics,'function_handle'));
verifyTrue(testCase, isa(problem.sim.functions,'function_handle'));

delete(cleanup);
end

function cleanup = addExamplePath(pathToAdd)
addpath(pathToAdd);
cleanup = onCleanup(@() rmpath(pathToAdd));
end
