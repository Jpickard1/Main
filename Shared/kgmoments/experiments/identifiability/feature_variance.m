%% Kronecker feature variances
% In this experiment, we look at the variance of different Kronecker
% features.

close all;
set(0,'DefaultLineLineWidth',1.1);
set(0,'DefaultAxesLineWidth',0.7);
set(0,'DefaultAxesFontSize',12);

%% Load experimental results
load coinflip_identify;

%%
nhist = 5; % number of histogram bins
Ti = 4;
trueF = zeros(nrep,4);
estF = zeros(nrep,4);
fitF = zeros(nrep,4);

for i=1:length(results)
    if results(i).Ti == Ti
        T = results(i).T; % assumes fixed T
        k = results(i).k; % assumes fixed k
        trueF(results(i).rep,1) = results(i).props.nedges;
        trueF(results(i).rep,2) = results(i).props.nwedges;
        trueF(results(i).rep,3) = results(i).props.ntripins;
        trueF(results(i).rep,4) = results(i).props.ntris;
        estF(results(i).rep,1) = results(i).stats.fitdata.nedges;
        estF(results(i).rep,2) = results(i).stats.fitdata.nwedges;
        estF(results(i).rep,3) = results(i).stats.fitdata.ntripins;
        estF(results(i).rep,4) = results(i).stats.fitdata.ntris;
        fitF(results(i).rep,1) = results(i).stats.sample_gdata.nedges;
        fitF(results(i).rep,2) = results(i).stats.sample_gdata.nwedges;
        fitF(results(i).rep,3) = results(i).stats.sample_gdata.ntripins;
        fitF(results(i).rep,4) = results(i).stats.sample_gdata.ntris;
        name = results(i).props.name;
    end
end

e = expected_kronecker_moments(T,k);
F(1) = e.nedges;
F(2) = e.nwedges;
F(3) = e.ntripins;
F(4) = e.ntris;

clf;
labels = {'edges','wedges','tripins','triangles'};
for i=1:4
    subplot(2,2,i);
    
    cla; hold on;
    [n,x] = hist(trueF(:,i),nhist);
    plot(x,n,'ko-');
    [n,x] = hist(estF(:,i),nhist);
    plot(x,n,'b.--');
    [n,x] = hist(fitF(:,i),nhist);
    plot(x,n,'rx-.');
    
    xl = xlim;
    %xx = linspace(xl(1),xl(2), 200);
    x1 = floor(xl(1)); x2 = ceil(xl(2)); xx = x1:ceil((x2-x1)/200):x2;
    %yy = normpdf(xx,F(i),sqrt(F(i)))/normpdf(F(i),F(i),sqrt(F(i)))*20;
    %yy = 20*poisspdf(xx,F(i))/poisspdf(round(F(i)),F(i));
    
    %plot(xx,yy,'-');
    
    ylim([0,25]);
    title(labels{i});
    
end
    

set_figure_size([4,4]);
print(sprintf('%s-features.eps',name),'-depsc2','-cmyk');