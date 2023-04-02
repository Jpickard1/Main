
maxN = 50;
xp = [];
poly = [];
for n=1:maxN
    xp = [xp; 2^n];
    poly = [poly; n^2];
end

x = 1:maxN;

figure('Renderer', 'painters', 'Position', [10 10 600 200]); hold on;
plot(x, xp);
plot(x, poly);
set(gca, 'YScale', 'log');
xlabel('Number of Vertices');
ylabel('Number of (Hyper)edges');
legend(["Hyperedges", "Edges"]);
saveas(gcf, 'complexity.png')

