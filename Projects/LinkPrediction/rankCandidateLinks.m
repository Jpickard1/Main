function [fEvalS, idx] = rankCandidateLinks(V, E, f, Ec)
%RANKCANDIDATELINKS This method ranks the likelihood of possible links
%forming in a network.
%
% PARAMS:
%       V: vertex set
%       E: edge set
%       f: link prediction function
%       Ec: set of candidate edges
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 9, 2023

% Construct adjacency matrix
A = real(E2A(V, E));

% Evaluate all candidate edges with f
fEval = zeros(size(Ec,1), 1);
for i=1:size(Ec,1)
    fEval(i) = f(A, Ec(i,1), Ec(i,1));
end

% Sort f
[fEvalS, idx] = sort(fEval, 'descend');

end
