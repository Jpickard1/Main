function A = small_world_alpha_model(n, k, alpha)
%SMALL_WORLD_ALPHA_MODEL This function constructs a small world graph
%   according to the procedure outlined in Chapter 3.1.1 (page 46-47) of
%   Small Worlds by Watts. The algorithms is designed to construct small
%   world graphs similar to how social networks are formed, where the
%   liklihood of 2 vertices being adjacent is a function of the number of
%   mutual neighbors.
%
%   NOTES:
%       - The algorithm terminated when k*n/2 edges have been chosen so
%       that the average degree is k.
%       - The book specifies that nodes 'choosing' the neighbor should be
%       selected in a random order such that no node chooses twice before
%       they all choose onces; however, for ease of implementation, because
%       the vertices are unlabeled, and because this only determines the
%       order in the last round through every vertex before the algorithm
%       terminates, this step is skipped and the nodes are selected in
%       normal order
%
% Auth: Joshua Pickard
% Date: May 11, 2022

%%
A = zeros(n,n);
terminate = false;

p = 0.1 * (nchoosek(n,2)^-1);
edges = 0;

while ~terminate
    for vx_i=1:n
        % Check for terminations
        if edges == k * n / 2
            terminate = true;
            break
        end
        neighbors_i = A(vx_i,:);
        % Compute R values
        R = zeros(n,1);
        for vx_j=1:n
            neighbors_j = A(vx_j,:);
            m = sum(neighbors_i .* neighbors_j);
            if m >= k
                r = 1;
            elseif m > 0
                r = (m/k)^alpha * (1-p) + p;
            else
                r = p;
            end
            R(vx_j) = r;
        end
        R(vx_i) = 0;
        P = R / sum(R);
        vx_chosen = randsample(length(P), 1, true, P);
        % Add in the edge
        A(vx_i, vx_chosen) = 1;
        A(vx_chosen, vx_i) = 1;
        edges = edges + 1;
    end
end

end

