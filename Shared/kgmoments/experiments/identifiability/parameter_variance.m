%% Kronecker feature variances
% In this experiment, we look at the variance of different Kronecker
% features.

addpath('~/dev/mcode');

close all;
set(0,'DefaultLineLineWidth',1.1);
set(0,'DefaultAxesLineWidth',0.7);
set(0,'DefaultAxesFontSize',12);

%% Load experimental results
load coinflip_identify;

%%
nhist = 6; % number of histogram bins
Ti = 4;
abcs = zeros(nrep,3);
bs = zeros(nrep,1);
cs = zeros(nrep,1);

maxr = 0;

for i=1:length(results)
    if results(i).Ti == Ti
        T = results(i).T; % assumes fixed T
        k = results(i).k; % assumes fixed k
        ri = results(i).rep;
        maxr = max(ri,maxr);
        abcs(ri,1) = results(i).params.a;
        abcs(ri,2) = results(i).params.b;
        abcs(ri,3) = results(i).params.c;
        name = results(i).props.name;
    end
end

%abcs = abcs(1:maxr,:);

clf;
labels = {'a','b','c'};
for i=1:3
    subplot(1,3,i);
    
    cla; hold on;
    % make a symmetric set
    x1 = T(i):0.005:T(i)+0.025;
    x = [fliplr(T(i)+(T(i)-x1(2:end))) x1];
    x(x>1) = [];
    x(x<0) = [];
    %[n,x] = hist(abcs(:,i),nhist);
    n = hist(abcs(:,i),x);
    plot(x,n,'rx-');
    stem(T(i),max(n));
    h=text(T(i),max(n)*(1.02+(1-max(n)/50)*0.13),sprintf('%s = %5.3f',labels{i},T(i)));
    set(h,'VerticalAlignment','bottom','HorizontalAlignment','center','Color','b');

%     xl = xlim;
%     %xx = linspace(xl(1),xl(2), 200);
%     x1 = floor(xl(1)); x2 = ceil(xl(2)); xx = x1:ceil((x2-x1)/200):x2;
%     %yy = normpdf(xx,F(i),sqrt(F(i)))/normpdf(F(i),F(i),sqrt(F(i)))*20;
%     yy = 20*poisspdf(xx,F(i))/poisspdf(round(F(i)),F(i));
%     
%     plot(xx,yy,'-');
    
    ylim([0,50]);
    if Ti==1 && i==3, xlim([0.2,0.3]); end;
    if Ti==2 && i==3, xlim([0.05,.15]); end;
    if Ti==2 && i==2, xlim([0.6,0.7]); end;
    if Ti==3 && i==2, xlim([0.2,0.3]); end;
    %title(labels{i});
end
    
%%


%%
set_figure_size([6,2]);
print(sprintf('%s-ident.eps',name),'-depsc2','-cmyk');