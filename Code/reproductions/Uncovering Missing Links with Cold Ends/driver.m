%% load matrix file
load('C:\Users\picka\Documents\my_projects\DBTM\MissingData\Data\USAir\adjacency_matrix.mat');
A = (A > 0.000);

%% Split the data
[training, testing] = partition_edges(A, 'rand', 0.8);

%% Compute similarity between each vertex
similarity = local_similarity(training, 'AA');

%% 
