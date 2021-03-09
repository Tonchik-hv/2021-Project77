function [idxSegm, structT, txtMsg] = calcSegms(varargin)
%--------------------------------------------------------
% the function calSegms(ts, 'parameter', parValue, ...) splits time series into periods
% Input:
% ts ([1, m]) - vector of observations, the time series to process 
% Parameters:
% 'T' - zero-order approximation of the average period-length
% 'start', 'end' - the indices specify the piece of the time series to be
% processed
% Output:
% idxSegm = matrix of logicals with the same number of columns as ts
%           0 means the point is inside the segment
%           1 means the point is the segment end
% structT - structure with fields 'mseT', 'fftT'.. contains the values of
% period length obtained by different methods
% txtTsLength - a message to output


txtMsg = [];
fignames = [];
ts = varargin{1};
if size(ts,1) < size(ts, 2)
    ts = ts';
end
T = [];
stp = 1;
endp = length(ts);
i = 2;
while i < nargin
    if strcmp(varargin{i}, 'T');
        T = varargin{i+1};
    elseif strcmp(varargin{i}, 'start');
        stp = varargin{i+1};
    elseif strcmp(varargin{i}, 'end');
        endp = varargin{i+1};
    elseif strcmp(varargin{i}, 'figname');
        fignames = varargin{i+1}; 
    end
    i = i + 2;
end    
ts = ts(stp:endp);
fftT = evalPeriod(ts, 'fft');
mseT =  evalPeriod(ts, 'mse');


if isempty(T)
    T = round(fftT);
end
%---- check if TS is long enought 
if 2*T >= endp - stp + 1
    endp = min(stp + 3*T + 1, length(ts));
    if 2*T >= endp - stp + 1 || ~isfinite(T)
        txtMsg = 'The time series are too short. The time series must contain at least three periods.';
        idxSegm = ~ismember(1:endp - stp + 1, [1:T:endp - stp + 1]);
        structT.mseT = mseT;
        structT.fftT = fftT;
        structT.segmT = T;
        structT.pcT = T;
        y = [0, 0];
        return;
    end
elseif 10*T < endp - stp + 1
    txtMsg = 'The time series are too long for demonstration. The piece of approximetely 10*T will be processed';
    ts = ts(stp:stp+10*T);
end

[dim1, dim2] = size(ts);
if dim2 > dim1
    ts = ts';    
end

   
[PCs, ~, lstLambda] = PCAnalysis(ts, round(1.5*T));
[y, pcT] = findPCPair(PCs, lstLambda, T);
if y(1) ~= 0
    idxSegm = PCSegmentation(PCs, y, pcT); %PCAnalysis(ts(:,i), periods(1), [1 2]);
else
    txtMsg = 'No periodical pair with T > 2 was detected.';
    idxSegm = [1 zeros(1, length(PCs)-2) 1];
end
idxPoints = find(idxSegm);
segmT = (idxPoints(end)-idxPoints(1))/(length(idxPoints)-1);


structT.mseT = mseT;
structT.fftT = fftT;
structT.segmT = segmT;
structT.pcT = pcT;

if ~isempty(fignames)
   plotPCs(PCs(y(1),:), PCs(y(2),:), idxPoints, fignames{1});
   plotLam(lstLambda, y, fignames{2});
end
end

function plotPCs(x, y, idxP, fig_name)
h = figure;
hold on;
plot(x, y, 'b-', 'markersize',8, 'linewidth', 2) % 1:length(idxSegm), tsSegm, 'b-', 
plot(x(idxP), y(idxP), 'ro', 'markersize', 8, 'linewidth', 2)
hold off;
xlabel('$\mathbf{y}_j$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
ylabel('$\mathbf{y}_{j+1}$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
set(gca, 'FontSize', 16, 'FontName', 'Times')
axis tight;
saveas(h, fig_name, 'png')
close(h);

end

function plotLam(lstLam, y, fig_name)
h = figure;
hold on;
plot(sqrt(lstLam), 'b-', 'markersize',8, 'linewidth', 2) % 1:length(idxSegm), tsSegm, 'b-',
if y(1)~= 0
plot(y, sqrt(lstLam(y)), 'ro', 'markersize', 8, 'linewidth', 2)
end
hold off;
xlabel('Ordered number, $j$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
ylabel('$\sqrt{\lambda_j}$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
set(gca, 'FontSize', 16, 'FontName', 'Times')
axis tight;
saveas(h, fig_name, 'png')
close(h);

end