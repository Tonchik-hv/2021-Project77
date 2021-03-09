function ARlen = estimFreqcy(ts)

% filtered ts
X = FFTAnalysis(ts, length(ts), [], 2:length(ts),0);


idxSegm = [];
[nX, xBars] = hist(X);
i = 2;
if sum(X < xBars(3)) < sum(X > xBars(8))
    sign = 1;
else
    sign = -1;
    nX = nX(end:-1:1);
end

while nX(i+1) > sum(nX(1:i)) %& i < 4
    %idxSegm = X*sign < sign*xBars(i);  %| X > xBars(11-i);
    i = i + 1;
    %idxSegm = X > par*max(X);
end

idxSegm = X*sign <= sign*xBars(i);
idxSegm = idxSegm';
idxSt = (idxSegm - [0 idxSegm(1:end-1)]) == 1;
idxEnd = (idxSegm - [idxSegm(2:end) 0]) == 1;

idxSt = find(idxSt.*length(X));
idxEnd = find(idxEnd.*length(X));

ARlen = length(X)/length(idxSt);
    


%idxSegm = X > par*max(X);
%{
idxSt = (idxSegm - [0 idxSegm(1:end-1)]) == 1;
idxEnd = (idxSegm - [idxSegm(2:end) 0]) == 1;

idxSt = find(idxSt.*length(X));
idxEnd = find(idxEnd.*length(X));

fs = SpS*length(idxSt)/length(X);
%}
end