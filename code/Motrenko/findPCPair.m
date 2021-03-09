function [y, T] = findPCPair(PCs, Lambda, T)


if nargin < 3
    T = [];
end
nPCs = sum(Lambda > 0);

idxPairs = [1:nPCs-1; 2:nPCs];
corrPy = [];
idx1 = [];
idx2 = [];
%c = [];
%r = [];
nFFT = length(PCs); %128;

for p = idxPairs
    [~, Pyy1, idxMax1] = FFTAnalysis(PCs(p(1),:)', 'nFFT', nFFT);
    [~, Pyy2, idxMax2] = FFTAnalysis(PCs(p(2),:)', 'nFFT', nFFT);
    idx1(end+1) = idxMax1(1);
    idx2(end+1) = idxMax2(1);
    corrPy(end+1) = corr(Pyy1, Pyy2);%*mean([Lambda(p(1)),Lambda(p(2))])/(Lambda(p(1)) - Lambda(p(2)));
    %[c(end+1), r(end+1)] = circlePars(PCs(p(1),:), PCs(p(2),:));
        
end
nP = find(abs(idx1 - idx2) < nFFT/100);
%trsh = mean(corrPy(abs(idx1 - idx2) < nFFT/100));
[Tpc, idxY] = max(nFFT*ones(size(nP))./(idx1(nP) - 1));
idxY = find(nFFT*ones(size(nP))./(idx1(nP) - 1) == Tpc);
[~, idx] = max(corrPy(nP(idxY)));
idxY = idxY(idx);

if isempty(idxY)
    y = [0,0];
    T = nFFT;
else
    y = idxPairs(:,nP(idxY))';
end

end


