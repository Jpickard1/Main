% See: Link Prediction in Complex Networks, Lu and Zhou, Equation 22
function RWR = randomWalkWithRestart_index(adj, i, j)
    c = 0.05;
    adj = adj / sum(adj(:)); % Row normalize adj matrix
    e_i = zeros(length(adj), 1);
    e_i(i) = 1;
    e_j = zeros(length(adj), 1);
    e_j(j) = 1;
    q_i = (1-c) * inv(1-(c*adj')) *e_i;
    q_j = (1-c) * inv(1-(c*adj')) *e_j;
    RWR = q_i(j) + q_j(i);
end
