%% Prelim Scratch
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 16, 2023

A = [1 1 0;
     1 1 1;
     0 1 1];

itrs = 3;
Akk = cell(3,1);

for i=1:itrs
    if i == 1
        Ak = kron(A, A);
    else
        Ak = kron(Ak, A);
    end
    [x, y] = find(Ak == 1);
    Akk{i} = [x, y];
    maxX = max(x);
    maxY = max(y);
end

figure('Renderer', 'painters', 'Position', [10 10 910 310]); hold on;
for i=1:itrs
    XY = Akk{i}; x = XY(:,1); y = XY(:,2);
    f = maxX / max(x); x = x .* f; y = y .* f;
    shift = (1.2 * maxX * (i - 1)); disp(shift);
    x = x + shift;
    scatter(x, y, '.b'); set ( gca, 'xdir', 'reverse' ); % xlim([0 maxX]); ylim([0 maxY]);
    % If i ~= 1 draw small box
    if i ~= 1
        rectangle('Position',[shift, 0, 28, 28])
    end
    % If i ~= 3 draw big box
    if i ~= 3
        rectangle('Position',[shift, 0, maxX+3, maxX+3]);
    end
end
% Draw lines connecting boxes
shift1 = 0;
shift2 = (1.2 * maxX * (1));
shift3 = (1.2 * maxX * (2));

X = [shift2 maxX+3]
Y1 = [28, maxX+3]
Y2 = [0 0]
inbetween = [Y1, fliplr(Y2)]
x2 = [X, fliplr(X)]
fill(x2, inbetween, [.7 .7 .7]);

X = [shift3 shift2+maxX+3]
x2 = [X, fliplr(X)]
fill(x2, inbetween, [.7 .7 .7]);

ylim([-2 maxY+4]);
xlim([-2 280]);
axis off;
set(gcf,'color','w');

saveas(gca, 'kronMultiscale.png')
