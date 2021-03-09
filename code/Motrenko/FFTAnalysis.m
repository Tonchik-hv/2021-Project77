function [filteredTS, Pyy, idxMax] = FFTAnalysis(varargin)  %(ts, N, SpS, idxFs, ifLog, ifPlot)
%--------------------------------------------------------------------------
% Implement Fourier-analysis using Fast Fourier Transformation 
% ts  - sudied time series, a vector if observations
% Nmax - number of frequencies to plot
% SpS - numeric SpS (samples per second) sets number of samples per second
%       Inf  - unset SpS, will be replaced by 1
%       [] - is used as signal to look for the first peak in periodigram
% idxN  -  indices of peak frequences (starting from the minimum frequency)
% modes:
% 'allMax' - return indices of all peaks
% 'idxPeaks' - return indicies of peaks, specified by idxN (idxN are the 
%              numbers of peaks oredered from the lowest fqc to the highest)
% 'maxPeaks' - return Nmax of most powerful peaks 


TRSH = 0.25;
ts = varargin{1};
% default values:
idxN = 1:length(ts); %indices of peak frequences (starting from the minimum frequency)
Nmax = length(ts); % max number of peaks (starting from the most powerfull peak)
nFFT = length(ts); 
%maxPeakidx = length(ts);
SpS = 1; 
idxFs = 1;
ifLog = 0; 
ifPlot = 0;

i = 2;
while i < length(varargin)
    switch varargin{i}
        case 'nFFT'
            nFFT = varargin{i+1};
            i = i+2;
        case 'allMax' % default: returns all detected peaks
            i = i+1;
        case 'idxPeaks' % returns peaks whose order is specified by idxPeaks
            idxN = varargin{i+1}; 
            i = i+2;
        case 'maxPeaks' % returns maxPeaks of peaks woth the greatest amplitudes
            Nmax = varargin{i+1};
            i = i+2;
        case 'ifPlot'
            ifPlot = varargin{i+1};
            i = i+2;
        case 'Frqs'
            Fs = varargin{i+1};
            i = i+2;
        case 'SpS'
            SpS = varargin{i+1};
            i = i+2;
        case 'ifLog'
            ifLog = varargin{i+1};
            i = i+2;    
    end
end

Y = fft(ts, nFFT);
Pyy = Y.*conj(Y)/length(ts);
Pyy = Pyy(1:round(length(ts)/2));
[idxAll, MaxVals] = findMax(Pyy, TRSH);
if idxAll == 1
    Pyy(1) = 0;
    [idxAll, MaxVals] = findMax(Pyy, TRSH);
end
%findMax_20Jan(Pyy, 0.1);

if Nmax < length(ts)
    [~, idxSorted] = sort(MaxVals, 'descend');
    idxSorted = idxSorted(1:min(length(idxAll), Nmax));
    idxMax = idxAll(idxSorted);
else
    idxMax = idxAll;
end
idxN = idxN(1:min(length(idxN), length(idxMax)));
idxMax = idxMax(idxN);


if ifPlot 
plotFFT(Pyy, idxMax, SpS, ifLog, TRSH);
end

%
filtrY = zeros(length(idxMax),1);
filtrY(idxMax) = Y(idxMax);
X = ifft(filtrY);
filteredTS = abs(X);
end


function [idxMax, tsMax] = findMax(X, trs)
if size(X,1) > size(X,2)
    X = X';
end


% Loess smothing:
XX = smooth(1:length(X),X,0.05,'loess')';
idxSegm = XX - min(XX(isfinite(XX))) > trs*(max(XX(isfinite(XX))))-min(XX(isfinite(XX)));


idxSt = (idxSegm - [0 idxSegm(1:end-1)]) == 1;
idxEnd = (idxSegm - [idxSegm(2:end) 0]) == 1;

idxSt = find(idxSt.*length(X));
idxEnd = find(idxEnd.*length(X));

%idxMax = zeros(1, length(X));
for i = 1:length(idxSt)
    [tsMax(i), idxMax(i)] = max(X(idxSt(i):idxEnd(i)));
    idxMax(i) = idxMax(i) + idxSt(i) - 1;       
end


end



function plotFFT(Pyy, idxMax, SpS, ifLog, trs)

f = SpS/length(Pyy)*(0:length(Pyy)-1);
if ~ifLog
figure;
plot(f,Pyy,'b-', f(idxMax), Pyy(idxMax), 'ro', 'linewidth', 2);
hold on;
plot(f, repmat(trs*max(Pyy), length(f),1), 'r-', 'linewidth', 2);
hold off;
xlabel('Frequency (Hz)', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
ylabel('Amplitude', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
set(gca, 'FontSize', 16, 'FontName', 'Times')
axis tight;
else
figure;
semilogy(f,Pyy,'b-', f(idxMax), Pyy(idxMax), 'ro', 'linewidth', 2);
hold on;
plot(f, repmat(trs*max(Pyy), length(f),1), 'r-', 'linewidth', 2);
hold off;
%title('Power spectral density')
xlabel('Frequency (Hz)', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
ylabel('Amplitude', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
set(gca, 'FontSize', 16, 'FontName', 'Times')
axis tight;
end

end