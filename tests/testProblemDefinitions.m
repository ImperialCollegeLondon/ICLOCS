function tests = testProblemDefinitions
%TESTPROBLEMDEFINITIONS Solver-free regression checks for example metadata.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(rootDir,'src')));
testCase.TestData.RootDir = rootDir;
end

function testBangBangFixedMeshDefinition(testCase)
exampleDir = fullfile(testCase.TestData.RootDir,'exampleProblems','BangBang','Fixed Mesh');
cleanup = addPathForTest(exampleDir);

[problem,guess] = BangBang;

verifyEqual(testCase, problem.time.t0_min, 0);
verifyEqual(testCase, problem.time.t0_max, 0);
verifyEqual(testCase, problem.time.tf_min, 0);
verifyEqual(testCase, problem.time.tf_max, 35);
verifyEqual(testCase, guess.tf, 30);
verifyEqual(testCase, problem.states.x0, [0 0]);
verifyEqual(testCase, problem.states.xfl, [300 0]);
verifyEqual(testCase, problem.states.xfu, [300 0]);
verifyEqual(testCase, problem.inputs.ul, -2);
verifyEqual(testCase, problem.inputs.uu, 1);
verifyTrue(testCase, isa(problem.analyticDeriv.gradCost,'function_handle'));

delete(cleanup);
end

function testOrbitRaisingOriginalDefinition(testCase)
exampleDir = fullfile(testCase.TestData.RootDir,'exampleProblems','OrbitRaising','Original_Formulation');
cleanup = addPathForTest(exampleDir);

[problem,guess] = OrbitRaising;

verifyEqual(testCase, problem.time.t0_min, 0);
verifyEqual(testCase, problem.time.tf_min, 3.32);
verifyEqual(testCase, problem.time.tf_max, 3.32);
verifyEqual(testCase, guess.tf, 3.32);
verifyEqual(testCase, problem.states.x0, [1 0 0 1]);
verifyEqual(testCase, problem.inputs.ul, [-1 -1]);
verifyEqual(testCase, problem.inputs.uu, [1 1]);
verifyEqual(testCase, problem.constraints.ng_eq, 1);
verifyEqual(testCase, problem.constraints.bl, [0 0]);
verifyEqual(testCase, problem.constraints.bu, [0 0]);
verifyEqual(testCase, problem.data.T1, 0.1405);
verifyEqual(testCase, problem.data.md, 0.0749);

delete(cleanup);
end

function cleanup = addPathForTest(pathToAdd)
addpath(pathToAdd);
cleanup = onCleanup(@() rmpath(pathToAdd));
end
