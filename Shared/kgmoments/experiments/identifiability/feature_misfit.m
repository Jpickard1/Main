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
Ti = 2;
trueF = zeros(nrep,4);
estF = zeros(nrep,4);
fitF = zeros(nrep,4);

for i=1:length(results)
    if results(i).Ti == Ti
        T = results(i).T; % assumes fixed T
        k = results(i).k; % assumes fixed k
        rep = results(i).rep;
        trueF(rep,1) = results(i).props.nedges;
        trueF(rep,2) = results(i).props.nwedges;
        trueF(rep,3) = results(i).props.ntripins;
        trueF(rep,4) = results(i).props.ntris;
        estF(rep,1) = results(i).stats.fitdata.nedges;
        estF(rep,2) = results(i).stats.fitdata.nwedges;
        estF(rep,3) = results(i).stats.fitdata.ntripins;
        estF(rep,4) = results(i).stats.fitdata.ntris;
        fitF(rep,1) = results(i).stats.sample_gdata.nedges;
        fitF(rep,2) = results(i).stats.sample_gdata.nwedges;
        fitF(rep,3) = results(i).stats.sample_gdata.ntripins;
        fitF(rep,4) = results(i).stats.sample_gdata.ntris;
        name = results(i).props.name;
    end
end
diffest = trueF - estF;
difffit = trueF - fitF;


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
    [n,x] = hist(diffest(:,i)./F(i),nhist);
    plot(x,n,'ko-');
    [n,x] = hist(difffit(:,i)./F(i),nhist);
    plot(x,n,'b.--');
    
    ylim([0,25]);
    title(labels{i});
    
end
    

set_figure_size([4,4]);
print(sprintf('%s-misfit.eps',name),'-depsc2','-cmyk');