function plotTS(ts, xTick)
if nargin < 2
    xTick = [];
end

color = {'b-', 'r-', 'g-'};
steps = []; %[80, 190, 310, 420];  %, 65, 90 ,115, 145, 170, 190, 219];
[n1, n2] = size(ts);
nTS = min(n1, n2); 
figure;
for i = 1:nTS
plot(ts(:,i) + 3*(i-1), color{i}, 'linewidth', 2);
hold on;
end
for i = 1:length(steps)
line([steps(i), steps(i)], [min(min(ts)), max(max(ts))+6]);
%hold on;
end
xlabel('Time $t$, s', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
ylabel('Acceleration', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
if ~isempty(xTick)
    set(gca, 'XTick',[1:length(ts)/(length(xTick)-1):length(ts)]);
    set(gca,'XTickLabel',xTick);
end
set(gca, 'FontSize', 16, 'FontName', 'Times')
axis tight;
hold off;
end