%% Check degree distribution smoothing with noise

%% add paths
addpath('~/dev/mcode');

%%

%P = [0.5 0.3 0.3];
%P = [0.9 0.5 0.5];
P = [0.57 0.19; 0.19 0.05]; % graph 500

k = 16;
ntrial = 100;

%%
% Show the problem in the degree distribution
figure(1); clf; hold on;
for t=1:ntrial
    a = P(1); b = P(2); c = P(3);
    G = rmat(k,[a b; b c],16);
    G = G|G';
    deg = sum(G,2);
    [n,x] = dhist(deg);
    loglog(n,x,'.'); 
    set(gca,'XScale','log'); set(gca,'YScale','log'); drawnow;
end

%%
% Show the "fix"
figure(1); clf; hold on;
df = zeros(2^k,1);
for t=1:ntrial
    a = P(1); b = P(2); c = P(3);
    G = rmatn(k,[a b; b c],1,16);
    G = G|G';
    deg = sum(G,2);
    [n,x] = dhist(deg);
    df(n(n>1)) = df(n(n>1)) + x(n>1);
    loglog(n,x,'.'); 
    set(gca,'XScale','log'); set(gca,'YScale','log'); drawnow;
end
loglog(df/ntrial,'r.');