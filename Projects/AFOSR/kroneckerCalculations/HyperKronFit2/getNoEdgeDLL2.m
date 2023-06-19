function DLL = getNoEdgeDLL2(theta, count, i, j, eLL)
%GETNOEDGEDLL
%
% Note: This code was generated using ChatGPT from a C++ function in the
% SNAP library and later modified by Joshua to fit into HAT.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: June 9, 2023

pDLL = count(i,j) / theta(i,j);
DLL = -theta(i,j) * pDLL / (1 - exp(eLL));

% DLL = (count(i,j) * (theta(i,j) ^ (count(i,j) - 1))) * (-1 * exp(eLL));

% DLL = -1 * count(i,j) * theta(i,j) / (1 - theta(i,j));
% DLL = -1 * count(i,j) * (theta(i,j) ^ (count(i,j) - 1)) / (1 - (theta(i,j) ^ count(i,j)));
% DLL = -1 * count(i,j) * (exp(eLL) / (theta(i,j))) / (1 - exp(eLL));

% DLL = sum(sum(theta .* count)) - theta(i,j);
% DLL = - count(i,j) * exp(DLL) / (1-exp(eLL));

end

% MtxDim = size(theta, 1);
% ThetaX = i; % mod(ParamId, GetDim());
% ThetaY = j; % fix(ParamId / GetDim());
% ThetaCnt = 0;
% DLL = 0;
% LL = 0;
% 
% for level = 0:NKronIters-1
%     X = mod(NId1, MtxDim);
%     Y = mod(NId2, MtxDim);
%     LVal = theta(X+1, Y+1);  % Assuming 1-based indexing in MATLAB
%     
%     if X == ThetaX && Y == ThetaY
%         if ThetaCnt ~= 0
%             DLL = DLL + LVal;
%         end
%         ThetaCnt = ThetaCnt + 1;
%     else
%         DLL = DLL + LVal;
%     end
%     
%     LL = LL + LVal;
%     NId1 = fix(NId1 / MtxDim);
%     NId2 = fix(NId2 / MtxDim);
% end
% 
% DLL = -ThetaCnt * exp(DLL) / (1.0 - exp(LL));