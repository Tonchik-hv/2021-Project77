function plotPeriods(ts, idxSegm, xTick, fig_name)

if nargin < 3
    xTick = [];
end
idxPoints = find(idxSegm);
h = figure;
hold on;
plot(1:length(ts), ts, 'b-', 'markersize',8, 'linewidth', 2) % 1:length(idxSegm), tsSegm, 'b-', 
plot(idxPoints, ts(idxPoints), 'ro', 'markersize', 8, 'linewidth', 2)
hold off;

if ~isempty(xTick)
    set(gca, 'XTick',[1:length(ts)/(length(xTick)-1):length(ts)]);
    set(gca,'XTickLabel',xTick);
end

xlabel('Time index $i$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
ylabel('Time series x(i)', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
set(gca, 'FontSize', 16, 'FontName', 'Times')
axis tight;
saveas(h, fig_name, 'png')
close(h);
end


