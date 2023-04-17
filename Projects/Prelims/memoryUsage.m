figure('Renderer', 'painters', 'Position', [0 0 500 500]); hold on;
hold on;
n = 10;
for i=1:n
    for j=1:n
        plot([i, i], [0, j], 'k');
        plot([0, j], [i, i], 'k');
    end
end
for i=1:n
    
end
title('Memory Access');
xlabel('Vertices'); ylabel('Vertices')
set(gca,'xtick',[]); set(gca,'ytick',[]);



