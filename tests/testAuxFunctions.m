function tests = testAuxFunctions
%TESTAUXFUNCTIONS Unit tests for low-level ICLOCS helper functions.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(rootDir,'src')));
testCase.TestData.RootDir = rootDir;
end

function testLinspaceMatInterpolatesColumnVectors(testCase)
actual = linspaceMat([1;10],[3;14],3);
expected = [1;2;3;10;12;14];
verifyEqual(testCase, actual, expected, 'AbsTol', 10*eps);
end

function testLinspaceMatReturnsEmptyForEmptyInput(testCase)
verifyEmpty(testCase, linspaceMat([],[],5));
end

function testNormalizeTMapsIntervalToMinusOneOne(testCase)
tau = normalizeT([0;5;10],0,10);
verifyEqual(testCase, tau, [-1;0;1], 'AbsTol', 10*eps);
end

function testSpkroneyeMatchesKronWithSparseIdentity(testCase)
A = sparse([1 0 2; 0 3 0]);
actual = spkroneye(3,A);
expected = kron(speye(3),A);
verifyEqual(testCase, full(actual), full(expected), 'AbsTol', 10*eps);
verifyTrue(testCase, issparse(actual));
end

function testLegendreEvalUsesLegendreBasis(testCase)
x = [-1;0;1];
coefficients = [1;2;0.5];
expected = 1 + 2*x + 0.5*((3*x.^2 - 1)/2);
actual = legendreEval(coefficients,x);
verifyEqual(testCase, actual, expected, 'AbsTol', 100*eps);
end

function testLegendrefitRoundTripsQuadraticData(testCase)
x = linspace(-1,1,9)';
y = 1 + 2*x + 0.5*((3*x.^2 - 1)/2);
coefficients = legendrefit(x,y,2,'qr');
actual = legendreEval(coefficients,x);
verifyEqual(testCase, actual, y, 'AbsTol', 1e-10);
end
