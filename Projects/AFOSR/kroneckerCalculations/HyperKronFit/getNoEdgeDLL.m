function DLL = getNoEdgeDLL(i, j, theta, NId1, NId2, NKronIters)
%GETNOEDGEDLL
%
% Note: This code was generated using ChatGPT from a C++ function in the
% SNAP library and later modified by Joshua to fit into HAT.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 9, 2023
    MtxDim = size(theta, 1);
    ThetaX = i; % mod(ParamId, GetDim());
    ThetaY = j; % fix(ParamId / GetDim());
    ThetaCnt = 0;
    DLL = 0;
    LL = 0;
    
    for level = 0:NKronIters-1
        X = mod(NId1, MtxDim);
        Y = mod(NId2, MtxDim);
        LVal = theta(X+1, Y+1);  % Assuming 1-based indexing in MATLAB
        
        if X == ThetaX && Y == ThetaY
            if ThetaCnt ~= 0
                DLL = DLL + LVal;
            end
            ThetaCnt = ThetaCnt + 1;
        else
            DLL = DLL + LVal;
        end
        
        LL = LL + LVal;
        NId1 = fix(NId1 / MtxDim);
        NId2 = fix(NId2 / MtxDim);
    end
    
    DLL = -ThetaCnt * exp(DLL) / (1.0 - exp(LL));
end
