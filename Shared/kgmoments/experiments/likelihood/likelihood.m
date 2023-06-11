addpath('../../matlab');
addpath('../../snap');
addpath('~/dev/matlab');

%%
A = readSMAT('../../data/ca-GrQc.smat');
A = A | A';
A = A - diag(diag(A));

%%
% Jure's parameters
%kronloglike(A,[0.999 0.245; 0.245 0.691])
kronloglike(A,[0.999 0.691; 0.681 0.245],1e6,1e6);

%%
% KronFit (rerun) parameters
kronloglike(A,[0.9999 0.3607; 0.3609 0.453],1e6,1e6);
kronloglike(A,[0.9999 0.350370; 0.350293 0.458220],1e6,1e6);

%%
% Our parameters
kronloglike(A,[1.000 0.467; 0.467 0.279],1e6,1e6)

%%
A = readSMAT('../../data/as-22july06.smat');
assert(nnz(diag(A)) == 0);
assert(isequal(A,A'));

%%
% Jure's parameters [0.954,0.594,0.019];
% Jure's log-likelihood (paper) -593,000 (ish)
kronloglike(A,[0.954 0.594; 0.594 0.019],1e6,1e6);

%%
% Our parameters
kronloglike(A,[1 0.611; 0.611 0],1e6,1e6);



% leskovec_kronfit_results('ca-HepPh') = [0.999,0.437,0.484];
% leskovec_kronfit_results('ca-HepTh') = [0.999,0.271,0.587];
% leskovec_kronfit_results('as20000102') = [0.987,0.571,0.049];
% leskovec_kronfit_results('as-22july06') = [0.954,0.594,0.019];