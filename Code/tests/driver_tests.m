function driver_tests()
%% DRIVER TESTS
%
% This script is the main script for testing, validating and updating the
% software for the Hardwired Genome project. This file should be run prior
% to any pull request or as frequently as needed when modifying the code.
%
% WRITTING NEW TESTS
%   Tests should be written to validate each level of the code and
%   organized as specified below. No level of test should call another
%   level of tests (i.e. system tests should not call the unit tests of the
%   associated code). Each test should be written as a function designed to
%   be called from this file. The types of tests, their scope, file and
%   data management are outlined below.
%
% TESTING FILE MANAGEMENT:
%   - Unit tests: unit tests should be written in the same directory as the
%   functions they test, saved with the following name test_unit_*....m
%   - System tests: system tests should test all of the components of a
%   subdirectory of Code/ and should be saved in this directory as
%   test_system_<directory name>.m
%   - Integration tests: Integration tests should test the usage of
%   multiple parts of the code together, and should be saved to this
%   directory as test_integration_<test name>.m
%   - Experimen tests: Experiment tests check multiple parts of the code
%   together similar to an integration test, but have the added purpose of
%   validating experimental results. These tests should be relavitly small
%   experiments, but they validate that the code produces precise (and
%   hopefully accurate) restuls.
%
% DATA:
%   There are 2 types of data that will be used in the above types of test.
%   1. Real data: real data will be used, saved, and handled similar to all
%      other data used in the repository
%   2. Test data: test data is sample, made up data for the purpose of
%   having checkable, smaller test cases we can run with handchecked
%   correct outputs. It will be saved in a similar file format in the
%   directory Data/tests
%
% Auth: Joshua Pickard (jpic@umich.edu)
% Date: May 23, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add paths
if ispc % DBTM office computer
    addpath(genpath("../utils"))
else % Great Lakes
    disp("Great Lakes needs to be set up")
end
add_paths()

% Store tests and results of tests
passed_tests = {};
failed_tests = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Termination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Testing summary:")

disp("Job Complete")
