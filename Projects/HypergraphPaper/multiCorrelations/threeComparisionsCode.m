function createfigure(X1, YMatrix1, X2, YMatrix2, YMatrix3)
%CREATEFIGURE(X1, YMatrix1, X2, YMatrix2, YMatrix3)
%  X1:  vector of plot x data
%  YMATRIX1:  matrix of plot y data
%  X2:  vector of plot x data
%  YMATRIX2:  matrix of plot y data
%  YMATRIX3:  matrix of plot y data

%  Auto-generated by MATLAB on 01-Mar-2023 15:09:34

% Create figure
figure1 = figure;

% Create subplot
subplot1 = subplot(1,3,1,'Parent',figure1);
hold(subplot1,'on');

% Create multiple line objects using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',subplot1,'Marker','x');
set(plot1(1),'DisplayName','Drezner');
set(plot1(2),'DisplayName','Wang and Zheng');
set(plot1(3),'DisplayName','Taylor');

% Create ylabel
ylabel('Average Vertex Degree');

% Create xlabel
xlabel('Multi-correlation Threshold');

% Create title
title('Degree Distribution');

hold(subplot1,'off');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,'Location','southwest');

% Create subplot
subplot2 = subplot(1,3,2,'Parent',figure1);
hold(subplot2,'on');

% Create multiple line objects using matrix input to plot
plot(X2,YMatrix2,'Parent',subplot2,'Marker','x');

% Create ylabel
ylabel('Number of Connected Components');

% Create xlabel
xlabel('Multi-correlation Threshold');

% Create title
title('Hypergraph Connectivity');

hold(subplot2,'off');
% Create subplot
subplot3 = subplot(1,3,3,'Parent',figure1);
hold(subplot3,'on');

% Create multiple line objects using matrix input to plot
plot(X1,YMatrix3,'Parent',subplot3,'Marker','x');

% Create ylabel
ylabel('Tensor Entropy');

% Create xlabel
xlabel('Multi-correlation Threshold');

% Create title
title('Entropy');

hold(subplot3,'off');