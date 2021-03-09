function idx = findMaxHist(X)

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

idx = X*sign <= sign*xBars(i);

end