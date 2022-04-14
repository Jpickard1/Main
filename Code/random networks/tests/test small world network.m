%% TEST SMALL WORLD NETWORK
%   This file contains a few basic test cases to ensure that the small
%   world networks are being generated correctly. A small world network is
%   parameterized with V, k (the mean degree), and beta (the rewiring
%   parameter). The initial version of 'small_world_network.m' created
%   these networks according the the Watts-Strogatz-Beta algorithm. I
%   modified thhis algorithm so that k should be a non-integer value, as
%   this allows for clearer comparisons between these and other types of
%   random networks.
%
%   4/14/2022 - The current tests included in this file are simple tests
%   where a random seed is set and the networks are generated. The
%   generated network is compared to known to be correct networks. This
%   test both the SW function as well as the wrapper function for
%   generating all kinds of networks. Joshua Pickard...
%
%   TODO: Finish adding some test cases
%         fix how edge rewiring works

%% 4/14/2022
params = containers.Map;
params('Beta') = 0;
V = 10;
p = 0.5;
network = generate_random_network('SW', V, p, params);
