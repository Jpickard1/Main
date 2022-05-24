function shotgun(A)
%SHOTGUN Call this function first for any data matrix. It produces a few
%   plots that are helpful for immediatly understanding data.
%
%   PLOTS
%   1. scree plot
%   2. degree distribution
%   3. graph (according to hall's drawing algorithm)
%   4. Fiedler vector
%
% Auth: Joshua Pickard jpic@umich.edu
% Date: May 23, 2022

figure;

% scree plot
subplot(2,2,1); scree(A);

% degree distribution
subplot(2,2,2); plot(sort(sum(A), 'descend'));
title('Degree Distribution'); xlabel('Vertex'); ylabel('Degree');

% Halls graph drawing algorithm
[v,~] = eigs(A, 2, 'smallestabs');
v = real(v);
subplot(2,2,3); gplot(A, v);
title('Graph');

% Fiedler vector
y1 = zeros(1,length(v))';
y2 = v(:,2);
x=1:length(y1);
x = [x;x];        % repeat x values
yy = [y1,y2]';   % vector of upper & lower boundaries
subplot(2,2,4); fill(x,yy,'b')    % fill area defined by x & yy in blue
title('Fiedler Vector');

end