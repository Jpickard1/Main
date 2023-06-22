clear; clc; close all;
path2data = "C:\Joshua\MissingData\Projects\KroneckerHypergraphPaper\3CAnalysis\";
simulationParameters = ["OFF", "ON", "HIGH"];
T = cell(200, 3);
missing = 0;
for s=1:3
    for i=1:200
        dataFile = path2data + "HG2theta/HKF_" + simulationParameters(s) + "_" + string(i) + "_0.01.mat";
        if isfile(dataFile)
            load(dataFile);
            T{i,s} = theta;
            disp(min(likelihoods));
            disp(max(likelihoods));
        else
            missing = missing + 1;
        end
    end
end
disp(missing);

E = cell(size(T));
for s=1:3
    for i=1:200
        disp(i);
        if ~isempty(T{i,s})
            theta = T{i,s};
            E{i,s} = heig(theta);
        end
    end
end
F = [];
O = [];
H = [];
for i=1:200
    if ~isempty(E{i,1}); F = [F std(E{i,1})]; end;
    if ~isempty(E{i,2}); O = [O std(E{i,2})]; end;
    if ~isempty(E{i,3}); H = [H std(E{i,3})]; end;
end

figure; hold on;
plot(sort(F));
plot(sort(O));
plot(sort(H));

%% Try norms
N = cell(size(T));
for s=1:3
    for i=1:200
        disp(i);
        if ~isempty(T{i,s})
            theta = T{i,s};
            N{i,s} = sum(theta,'all');
        end
    end
end
F = [];
O = [];
H = [];
for i=1:200
    if ~isempty(N{i,1}); F = [F (N{i,1})]; end;
    if ~isempty(N{i,2}); O = [O (N{i,2})]; end;
    if ~isempty(N{i,3}); H = [H (N{i,3})]; end;
end

figure; hold on;
plot(sort(F));
plot(sort(O));
plot(sort(H));
