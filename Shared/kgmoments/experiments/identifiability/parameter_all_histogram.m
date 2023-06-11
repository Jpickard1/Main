% I decided not to use this histogram as it compresses the information too
% much to be useful.  The individual plots can show the detail around the
% fit better.

clf; cla; hold on;
labels = {'a','b','c'};
linstrs = {'k.-','bx-','r*-'}
for i=1:3
    
    % make a symmetric set
    x1 = T(i):0.005:T(i)+0.02;
    x = [fliplr(T(i)+(T(i)-x1(2:end))) x1];
    x(x>1) = [];
    x(x<0) = [];
    %[n,x] = hist(abcs(:,i),nhist);
    n = hist(abcs(:,i),x);
    plot(x,n,linstrs{i});
    stem(T(i),max(n),[linstrs{i}(1) '-']);
    h=text(T(i),max(n)*1.1,labels{i});
    set(h,'VerticalAlignment','bottom','HorizontalAlignment','center');

%     xl = xlim;
%     %xx = linspace(xl(1),xl(2), 200);
%     x1 = floor(xl(1)); x2 = ceil(xl(2)); xx = x1:ceil((x2-x1)/200):x2;
%     %yy = normpdf(xx,F(i),sqrt(F(i)))/normpdf(F(i),F(i),sqrt(F(i)))*20;
%     yy = 20*poisspdf(xx,F(i))/poisspdf(round(F(i)),F(i));
%     
%     plot(xx,yy,'-');
    
    ylim([0,40]);
    xlim([-0.05,1.05]);
 
end