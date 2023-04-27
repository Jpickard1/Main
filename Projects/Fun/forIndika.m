close all

theta0 = [0.01 2*pi/3 4*pi/3];
tspan = [0 8]; 

f = @(~,theta) sum(sin(theta'-theta),2);

[t, theta] = ode45(f, tspan, theta0);

figure(1)
plot(t, theta,'LineWidth',1.5);
xlabel("Time"), ylabel("\theta")
legend(["\theta_1","\theta_2","\theta_3"])
set(gca,'TickLabelInterpreter','latex',...
    'YTick',theta0,...
    'YTickLabel',{'0','$\frac{2\pi}{3}$','$\frac{4\pi}{3}$'});
axis square

figure(2)
plot3(theta(:, 1), theta(:, 2), theta(:, 3),'LineWidth',1.5);
xlabel("\theta_1"), ylabel("\theta_2"), ylabel("\theta_3")
axis square
grid on


figure(3)
n = 10;
[Theta1, Theta2, Theta3] = meshgrid(0:pi/n:5*pi/2);
dTheta1 = sin(Theta2-Theta1)+sin(Theta3-Theta1);
dTheta2 = sin(Theta1-Theta2)+sin(Theta3-Theta2);
dTheta3 = sin(Theta1-Theta3)+sin(Theta2-Theta3);
quiver3(Theta1, Theta2, Theta3, dTheta1, dTheta2, dTheta3, 'r');

lim = [0 5*pi/2];
tix = [0 pi/2 pi 3*pi/2 2*pi 5*pi/2];
tixlabels = ...
    {'0','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$','$\frac{5\pi}{2}$'};

set(gca,'TickLabelInterpreter','latex', ...
    'XLim',lim,'XTick',tix,'XTickLabel',tixlabels, ...
    'YLim',lim,'YTick',tix,'YTickLabel',tixlabels, ...
    'ZLim',lim,'ZTick',tix,'ZTickLabel',tixlabels, ...
    'LineWidth',1.5)
axis square

%{
% Add stream line from Theta1=0, Theta2=pi, Theta3=pi
hold on
streamline(Theta1, Theta2, Theta3, dTheta1, dTheta2, dTheta3, 0, pi, pi);
s = gca;
s.Children(1).LineWidth = 1.5;

% Find equally spaced time points

numTimePoints = 3;
theta1Vals = linspace(0, 1.5, numTimePoints);
idx = zeros(1, numTimePoints);

for i = 1:numTimePoints
    [~, idx(i)] = min(abs(theta(:, 1)-theta1Vals(i)));
end

text(theta(idx,1), theta(idx,2), theta(idx,3)-[0.5;0.3;0], arrayfun(@(x) sprintf("t=%.2f",x), tspan(idx)), 'FontSize', 12)
hold off
%}