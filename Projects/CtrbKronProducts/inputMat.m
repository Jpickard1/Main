function [B] = inputMat(ctrls, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B = zeros(n, numel(ctrls));
for i=1:numel(ctrls)
    B(ctrls(i),i) = 1;
end
end

