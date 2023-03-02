function symVec = loadSymVecs(n, l)
%LOADSYMVECS This function loads symbolic vectors from .mat files.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 26, 2023
    fName = "symVecs/" + string(n) + "_" + string(l) + string(".mat");
    load(fName)
end