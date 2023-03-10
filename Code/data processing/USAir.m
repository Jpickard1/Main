% Description: This file makes an adjacency matrix from a pairs file. A
% pairs file is a file with 3 columns: node i, node j, and the edge
% connecting them
%
% file_path = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Data\USAir\edge_pairs.txt'
% 
in = tdfread(file_path, '\t');
data = in.x_v_v_e;
v = max([data(:,1); data(:,2)]);
A = zeros(v);
for e=1:height(data)
    i = data(e,1);
    j = data(e,2);
    edge_weight = data(e,3);
    A(i,j) = A(i,j) + edge_weight;
    A(j,i) = A(j,i) + edge_weight;
end
save C:\Users\picka\Documents\my_projects\DBTM\Main\Data\USAir\adjacency_matrix.mat A


% Whenever you get data compute these 3 things:
% Fiddler number
% Condition number
% Entropy
