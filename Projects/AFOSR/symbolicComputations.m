%% Symbolic Computations
%
%   The purpose of this file is to precompute the symbolic vectors required
%   to select the minimal observable node sets. Constructing large symbolic
%   vectors is expensive, but their multiplication with sparse matrices is
%   relatively cheap. To improve the speed of analyzing many small
%   hypergraphs (i.e. n<=7) I precompute all possible required symbolic
%   vectors in this script.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 22, 2023
function symbolicComputations(n)
    disp(n)
    maxSize = 1e15;
    
    f = dir(fullfile('symVecs', '*.mat'));

    x = sym('x_%d',[n 1]);      % Set symbolic state vector

    % Get largest symbolic vector for n that has been precomputed
    i = 0;
    maxSizeP = 0;
    maxFile = "";
    for j=1:length(f)
        fName = f(j).name;
        % Check if this has been computed yet
        if startsWith(fName, string(n) + "_")
            symSize = str2num(fName(3:end-4));
            if symSize > maxSizeP
                i = log(symSize) / log(n);
                maxSizeP = symSize;
                maxFile = fName;
            end
        end
    end
    if i ~= 0
        disp("symVecs/" + string(maxFile));
        load("symVecs/" + string(maxFile));
    else
        symVec = sym('x_%d',[n 1]);      % Set symbolic state vector
    end
    disp("Loaded from File");
    disp(length(symVec));
	disp(maxSizeP);
	disp(maxSize);
	while length(symVec) < maxSize
        if i~= 0
            disp(i);
            symVec = kronSymVec(symVec, x, i-1);
        end
        disp(length(symVec));
        fileName = "symVecs/" + string(n) + "_" + string(length(symVec)) + ".mat";
        cmd = "save " + fileName + " " + "symVec -v7.3";
        eval(cmd);
        disp(cmd);
        disp(i); i = i + 1;
    end 
end 

%% kronSymVec

function symVec = kronSymVec(symVec, x, pow)
    disp("    kronSymVec call: " + string(length(x)) + "-" + string(pow));
    n = length(x);
    totalSize = n^pow;
    while length(symVec) < totalSize
        symVec = kron(symVec, x);
    end
end

