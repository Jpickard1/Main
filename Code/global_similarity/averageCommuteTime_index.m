% See: Link Prediction in Complex Networks, Lu and Zhou, Equation 19
function ACT = averageCommuteTime_index(adj, i, j)
    deg = adj2deg(adj);
    L = deg - adj;
    L_i = pinv(L_i);
    l_ii = L_i(i,i);
    l_jj = L_i(j,j);
    l_ij = L_i(i,j);
    ACT = 1 / (l_ii + l_jj + (2*l_ij));
end
