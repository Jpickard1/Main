function Output = DMD(DataMatrix,MetaData,thresh)

% Code by Joshua Proctor at the Institute for Disease Modeling 2015
% Please email joshlproctor@gmail.com for appropriate citations and questions


% disp('Performing DMD')
% disp('')
% disp('')

% Construct data matrices note, time snapshots are columns
    X  = DataMatrix(:,1:end-1);  
    Xp = DataMatrix(:,2:end);   

% Compute the SVD 
    [UX,SigX,VX] = svd(X,'econ');
    
% Choose the number of singular values based on energy percent
    r = find(cumsum((diag(SigX)/sum(diag(SigX)))) > thresh,1);
    
% Compute A_tilde    
    A = UX(:,1:r)'*Xp*VX(:,1:r)*inv(SigX(1:r,1:r));
    
% Compute Eigenanalysis
    [V,D] = eig(A);

% Compute Dynamic Modes
    Dyn = Xp*VX(:,1:r)*inv(SigX(1:r,1:r))*V;

% Place computed variables in a single output file 
Output.DataMatrix = DataMatrix;
Output.MetaData   = MetaData;

Output.X = X;
Output.Xp   = Xp;
     
Output.DMD.D = diag(D);
Output.DMD.V = V;
Output.DMD.DynamicModes = Dyn;
Output.DMD.A = A;
Output.DMD.Sig = SigX;
Output.DMD.r = r;
Output.DMD.UX = UX;
Output.DMD.VX = VX;         

%if size(DataMatrix,1) <= 1000
% Compute A_bar SCOTT ADD
    A_bar = Xp*VX*pinv(SigX)*UX';
    [V_bar,D_bar] = eig(A_bar);
    
    Output.DMD.A_bar = A_bar;           %SCOTT ADD
    Output.DMD.D_bar = diag(D_bar);     %SCOTT ADD
    Output.DMD.V_bar = V_bar;           %SCOTT ADD
%end
end