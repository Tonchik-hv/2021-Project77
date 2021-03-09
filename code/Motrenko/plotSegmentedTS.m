function plotSegmentedTS(ts, idxSegm, idxTrue, mode)

time = 1:length(ts);
if nargin < 4   
    idxClst = unique(idxSegm);
    
    color = {'b-', 'r-', 'g-','m-', 'c-', 'k-','y-'};
    nClst = 0;
    figure;
    for clst = idxClst'
        nClst = nClst + 1;
    tsSegm = ones(1,length(idxSegm))*NaN;
    tsSegm(idxSegm == clst) = ts(idxSegm == clst);

    %time(idxSegm == clst)
    plot(1:length(tsSegm), tsSegm, color{nClst}, 'linewidth', 2)
    hold on    
    end
    for i = 1:length(idxTrue)
    line([idxTrue(i), idxTrue(i)],[min(ts), max(ts)])
    hold on
    end
else 
    figure;
    plot(time, ts, 'b-', 'linewidth', 2)
    hold on;
    for i = 1:length(idxSegm)
        line([idxSegm(i), idxSegm(i)],[min(ts), max(ts)])
    end
end
xlabel('Time, $t$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
ylabel('Value', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
set(gca, 'FontSize', 16, 'FontName', 'Times')
axis tight;
hold off
end