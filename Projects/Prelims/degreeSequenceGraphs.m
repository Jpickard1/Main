%% DEGREE SEQUENCE GRAPHS
%
%   Here I generate plots to compare the number of possible graphs on n
%   labeled vertices under 3 conditions:
%
%   1. No information is known about the edges
%   2. The total number of edges is known
%   3. The degree sequence is fully known
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 12, 2023

maxN = 40;
x = [];
pct = 0:1/maxN:1;
N = 2:maxN;
C1 = zeros(length(pct), length(N));
C2 = zeros(length(pct), length(N));
C3 = zeros(length(pct), length(N));
for ni=1:length(N)
    n = N(ni);
    for mi=1:length(pct)
        m = pct(mi);
        M = round(m * nchoosek(n, 2));
        C1(mi, ni) = 2^nchoosek(n,2);
        C2(mi, ni) = nchoosek(nchoosek(n, 2), M);
        C3(mi, ni) = nchoosek(nchoosek(n, 2), M) / nchoosek(n + M, M);
    end
end

C3 = 100 * C3 ./ max(max(C1));
C2 = 100 * C2 ./ max(max(C1));
C1 = 100 * C1 ./ max(max(C1));

c1 = zeros(length(pct), length(N), 3);
c2 = zeros(length(pct), length(N), 3);
c3 = zeros(length(pct), length(N), 3);
for i=1:length(pct)
    for j=1:length(N)
        c1(i,j,:) = [1 0 0];
        c2(i,j,:) = [0 1 0];
        c3(i,j,:) = [0 0 1];
    end
end


X = repmat(N, length(pct), 1);
Y = repmat(pct', 1, length(N));
figure; hold on;
surf(X, Y, C1, c1);
surf(X, Y, C2, c2);
surf(X, Y, C3, c3);
%set(gca, 'ZScale', 'log');
xlabel('Number of Vertices');
ylabel('Number of Edges');
zlabel('Number of Graphs');
legend(["No Information", "Known Number of Edges", "Known Degree Sequence"])


%% 
figure; hold on; 
set(gca, 'YScale', 'log');
scatter(2:maxN, C1, '.');
scatter(2:maxN, C2, '.');
scatter(2:maxN, C3, '.');
% scatter(x, C1, '.');
% scatter(x, C2, '.');
% scatter(x, C3, '.');
title('Utility of Edge Distribution in Bounding Number of Graphs')
xlabel('Number of Vertices'); ylabel('Number of Possible Graphs')
legend(["No Information", "Known Number of Edges", "Known Degree Sequence"])

%% 
maxN = 40;
C1 = [];
C2 = [];
C3 = [];
x = [];
pct = 0.2;
for n=2:maxN
    c1 = 2^nchoosek(n,2);
    c2 = nchoosek(nchoosek(n,2), round(pct * nchoosek(n,2)));
    c3 = nchoosek(nchoosek(n,2), round(pct * nchoosek(n,2))) / nchoosek(n + round(pct * nchoosek(n,2)), round(pct * nchoosek(n,2)));
    C1 = [C1; c1];
    C2 = [C2; c2];
    C3 = [C3; c3];
    x = [x; nchoosek(n, 2)];
end
figure; hold on; 
set(gca, 'YScale', 'log');
scatter(2:maxN, C1, '.');
scatter(2:maxN, C2, '.');
scatter(2:maxN, C3, '.');
% scatter(x, C1, '.');
% scatter(x, C2, '.');
% scatter(x, C3, '.');
title('Utility of Edge Distribution in Bounding Number of Graphs')
xlabel('Number of Vertices'); ylabel('Number of Possible Graphs')
legend(["No Information", "Known Number of Edges", "Known Degree Sequence"])

%% Complexity Classes Plotted on a Log Scale

maxN = 40;
C1 = [];
C2 = [];
C3 = [];
pct = 0.2;
for n=2:maxN
    c1 = factorial(n);
    c2 = 5^n;
    c3 = n^5;
    C1 = [C1; c1];
    C2 = [C2; c2];
    C3 = [C3; c3];
end
figure; hold on; 
set(gca, 'YScale', 'log');
scatter(2:maxN, C1, '.');
scatter(2:maxN, C2, '.');
scatter(2:maxN, C3, '.');
title('Utility of Edge Distribution in Bounding Number of Graphs')
xlabel('Number of Vertices'); ylabel('Number of Possible Graphs')
legend(["Factorial", "Exponential", "Polynomial"])
