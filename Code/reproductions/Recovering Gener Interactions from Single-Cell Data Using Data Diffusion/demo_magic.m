%% D
num_cells = 20;
num_genes = 500;
D = rand(num_cells, num_genes);

known = randi([0 9], num_cells, num_genes);
known = (known > 0);
D_known = D;
D_known(~known) = 0;

%% Dist
dist = pdist(D);
dist = squareform(dist);
%% A 
sigma = 5;
A = zeros(size(dist));
for i=1:length(A)
    for j=1:length(A)
        A(i,j) = exp(-(dist(i,j)/sigma)^2);
    end
end
%% M
M = zeros(size(A));
for i=1:length(A)
    M(i,:) = A(i,:) / sum(A(i,:));
end

%% Impute the data
max_t = 30;
imputed_measured_norms = zeros(max_t, 1);
imputed_rank = zeros(max_t, 1);
for t=1:max_t
    if t ~= 1
        M_exp = M_exp * M;
    else
        M_exp = M;
    end
    D_imputed = M_exp * D;
    r = rank(D_imputed, 1e-4);
    imputed_rank(t) = r;
    % Mask the unknown values
    D_imputed_known = D_imputed;
    D_imputed_known(~known) = 0;
    imputed_measured_norm = norm(D_imputed_known - D_known, 'fro');
    imputed_measured_norms(t) = imputed_measured_norm;

end
figure;
subplot(2,1,1);
plot(1:max_t, imputed_rank);
xlabel('t');
ylabel('Imputed rank');
subplot(2,1,2);
plot(1:max_t, imputed_measured_norms);
xlabel('t');
ylabel('Norm');


