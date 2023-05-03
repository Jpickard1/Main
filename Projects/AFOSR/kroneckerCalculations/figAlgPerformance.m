%% Algorithm Performance
%
%   This file generates some plots to describe the performance of my NTKP
%   algorithm
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: May 2, 2023

%% Set experiment parameters
itrs = 1;
maxN = 15;

% Set recording parameters
newError = containers.Map();
oldError = containers.Map();
ref = containers.Map();

newError("Time") = []; % zeros(length(2:maxN), itrs); %[];
newError("N1")   = []; % zeros(length(2:maxN), itrs); %[];
newError("NE")   = []; % zeros(length(2:maxN), itrs); %[];
newError("NI")   = []; % zeros(length(2:maxN), itrs); %[];
newError("NP")   = []; % zeros(length(2:maxN), itrs); %[];
newError("SV")   = []; % zeros(length(2:maxN), itrs); %[];
oldError("Time") = []; % zeros(length(2:maxN), itrs); %[];
oldError("N1")   = []; % zeros(length(2:maxN), itrs); %[];
oldError("NE")   = []; % zeros(length(2:maxN), itrs); %[];
oldError("NI")   = []; % zeros(length(2:maxN), itrs); %[];
oldError("NP")   = []; % zeros(length(2:maxN), itrs); %[];
oldError("SV")   = []; % zeros(length(2:maxN), itrs); %[];
ref("N1")        = []; % zeros(length(2:maxN), itrs); %[];
ref("NE")        = []; % zeros(length(2:maxN), itrs); %[];
ref("NI")        = []; % zeros(length(2:maxN), itrs); %[];
ref("NP")        = []; % zeros(length(2:maxN), itrs); %[];
ref("SV")        = []; % zeros(length(2:maxN), itrs); %[];

%% experiment with random tensors
for n=2:maxN
    for i=1:itrs
        % Make tensor
        A = rand(n^2,n^2,n^2);
        
        % Perform and time 2 tensor decompositions
        tic;        
        [A1, A2] = NTKP(A);
        newError("Time") = [newError("Time") toc];
        tic; 
        [B, Sigma] = tkpsvd(A, n * ones(1,6), 1); 
        oldError("Time") = [oldError("Time") toc];

        % Reconstruct tensors
        NEW = superkron(A1, A2);        
        OLD = Sigma(1) * superkron(B{1,1}, B{2,1});
        
        % identify differences between new and old tensors
        N = A - NEW;
        O = A - OLD;

        % Compute tensor norms and save results
        [N1, NE, NI, NP] = tensorNorms(N);
        newError("N1") = [newError("N1") N1];
        newError("NE") = [newError("NE") NE];
        newError("NI") = [newError("NI") NI];
        newError("NP") = [newError("NP") NP];
        [N1, NE, NI, NP] = tensorNorms(O);
        oldError("N1") = [oldError("N1") N1];
        oldError("NE") = [oldError("NE") NE];
        oldError("NI") = [oldError("NI") NI];
        oldError("NP") = [oldError("NP") NP];
        
        % Save reference values
        [N1, NE, NI, NP] = tensorNorms(A);
        ref("N1") = [ref("N1") N1];
        ref("NE") = [ref("NE") NE];
        ref("NI") = [ref("NI") NI];
        ref("NP") = [ref("NP") NP];

        % Check singular values
        Amat = reshape(A, [size(A,1), numel(A)/size(A,1)]);
        Omat = reshape(O, [size(A,1), numel(A)/size(A,1)]);
        Nmat = reshape(N, [size(A,1), numel(A)/size(A,1)]);
        A1 = svds(Amat, 1);
        O1 = svds(Omat, 1);
        N1 = svds(Nmat, 1);
        oldError("SV") = [oldError("SV") abs(A1 - O1)];
        newError("SV") = [newError("SV") abs(A1 - O1)];
    end
    disp(n)
end

%% Experiment with random hypergraphs
for n=3:maxN
    for i=1:itrs
        % Make tensor
        e = nchoosek(n, 3); 0.25;
        H1 = HAT.uniformErdosRenyi(n, e, 3); A1 = H1.adjTensor;
        H2 = HAT.uniformErdosRenyi(n, e, 3); A2 = H2.adjTensor;
        A = superkron(A1, A2);
        
        % Perform and time 2 tensor decompositions
        tic;        
        [A1, A2] = NTKP(A);
        newError("Time") = [newError("Time") toc];
        tic; 
        [B, Sigma] = tkpsvd(A, n * ones(1,6), 1); 
        oldError("Time") = [oldError("Time") toc];

        % Reconstruct tensors
        NEW = superkron(A1, A2);        
        OLD = Sigma(1) * superkron(B{1,1}, B{2,1});
        
        % identify differences between new and old tensors
        N = A - NEW;
        O = A - OLD;

        % Compute tensor norms and save results
        [N1, NE, NI, NP] = tensorNorms(N);
        newError("N1") = [newError("N1") N1];
        newError("NE") = [newError("NE") NE];
        newError("NI") = [newError("NI") NI];
        newError("NP") = [newError("NP") NP];
        [N1, NE, NI, NP] = tensorNorms(O);
        oldError("N1") = [oldError("N1") N1];
        oldError("NE") = [oldError("NE") NE];
        oldError("NI") = [oldError("NI") NI];
        oldError("NP") = [oldError("NP") NP];
        
        % Save reference values
        [N1, NE, NI, NP] = tensorNorms(A);
        ref("N1") = [ref("N1") N1];
        ref("NE") = [ref("NE") NE];
        ref("NI") = [ref("NI") NI];
        ref("NP") = [ref("NP") NP];

        % Check singular values
        Amat = reshape(A, [size(A,1), numel(A)/size(A,1)]);
        Omat = reshape(O, [size(A,1), numel(A)/size(A,1)]);
        Nmat = reshape(N, [size(A,1), numel(A)/size(A,1)]);
        A1 = svds(Amat, 1);
        O1 = svds(Omat, 1);
        N1 = svds(Nmat, 1);
        oldError("SV") = [oldError("SV") abs(A1 - O1)];
        newError("SV") = [newError("SV") abs(A1 - O1)];
    end
    disp(n)
end

%% Plot results
x = (3:maxN) .^ 2;
nPlots = 6;
figure;
subplot(1, nPlots, 1); title('Run Time'); hold on;
plot(x, newError("Time"))
plot(x, oldError("Time"));
subplot(1, nPlots, 2,'YScale', 'log'); title('Norm 1'); hold on;
plot(x, newError("N1"))
plot(x, oldError("N1"));
plot(x, ref("N1"));
subplot(1, nPlots, 3); title('Euclidean Norm'); hold on;
plot(x, newError("NE"))
plot(x, oldError("NE"));
plot(x, ref("NE"));
subplot(1, nPlots, 4); title('Max Norm'); hold on;
plot(x, newError("NI"))
plot(x, oldError("NI"));
plot(x, ref("NI"));
subplot(1, nPlots, 5,'YScale', 'log'); title('Inner Product (Squared Euclidean)'); hold on;
plot(x, newError("NP"))
plot(x, oldError("NP"));
plot(x, ref("NP"));
legend(["New", "Old", "True"]);
subplot(1, nPlots, 6); title('Singular Values'); hold on;
plot(x, newError("SV"))
plot(x, oldError("SV"));
