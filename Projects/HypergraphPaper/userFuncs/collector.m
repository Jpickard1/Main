function [O] = collector(ds,varargin)
%COLLECT Summary of this function goes here
%   Detailed explanation goes here
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date; February 9, 2023

O = varargin;
filename = ds + "_Observations.mat";
save(filename, "O")

end

