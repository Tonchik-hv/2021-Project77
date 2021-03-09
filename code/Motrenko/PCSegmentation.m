function idxSegm = PCSegmentation(PCs, y, T, plotTs)

if nargin < 4 
    plotTs = [];
    
end
%oneCut(PCs(y(1),:),PCs(y(2),:), T);
idxSegm = newSegm(PCs(y(1),:),PCs(y(2),:), T);
%{
if length(PCs) > 5000
    phi = 0;
else
    phi = [];
end
%--- 18/12, prifiling
phi = 0;
[idxSegm, phi] = segment1(PCs(y,:), phi);
%}

%plotPCs(PCs(y,1:500),2, phi);
%idxSegm(idxMax) = [idxSegmM, zeros(1, sum(idxMax) - length(idxSegmM))];
%idxSegm = segment2(PCs(y,:));

if ~isempty(plotTs)
    plotPeriods(plotTs, idxSegm); % , [0:100:500]/20  [0:500:round(length(plotTs))]
end

end

function idxSegm = newSegm(x, y, T)

% Choose a good point for dissection:
K = convhull(x,y);
K = K(1:end-1);
if length(K) > T*0.5
     segmT = length(K);%sum(diff(K)<3);
    
    [~, k] = max(diff(K));
    idxClose = findClosePoints(x, y, k);
    while idxClose(end) < length(x) - segmT*1.1
        [idxClose, iNew] = fixClose(idxClose, length(idxClose),x, y, segmT);
        if iNew == 0
            break;
        end
    end
    while idxClose(1) > segmT*1.1
        [idxClose, iNew] = fixClose(idxClose, 0,x, y, segmT);
        if iNew == 0
            break;
        end
    end
    idxClose = idxClose(idxClose < length(x));
    segms = diff(idxClose);
    while max(segms)/median(segms) > 1.5
        [~, i] = max(segms);
        [idxClose, iNew] = fixClose(idxClose, i,x, y);
        segms = diff(idxClose);
    end
    if length(idxClose) < length(x)/T - 1
        idxClose = fixClose(idxClose, 0,x, y);
    end
else
var = [];
check = [];
for k = K'    
    %dist = sqrt((x-x(k)*ones(size(x))).^2 + (y-y(k)*ones(size(y))).^2);
    idxClose = findClosePoints(x, y, k);
    %segmT = median(diff(idxClose));
    while idxClose(end) < length(x) - T*1.1
        [idxClose, iNew] = fixClose(idxClose, length(idxClose),x, y, T);
        if iNew == 0
            break;
        end
    end
    while idxClose(1) > T*1.1
        [idxClose, iNew] = fixClose(idxClose, 0,x, y, T);
        if iNew == 0
            break;
        end
    end
    segms = diff(idxClose);
    while max(segms)/median(segms) > 1.5
        [~, i] = max(segms);
        [idxClose, iNew] = fixClose(idxClose, i,x, y);
        segms = diff(idxClose);
        if iNew == 0
            break;
        end
    end
    [var, check] = checkPoints(idxClose, length(x), T); 
    if isfinite(var) && check
        break;
    end
    %[var(end+1), check(end+1)]  = checkPoints(idxClose, length(x), T);        
end
end
%{
if sum(check) > 0 && sum(isfinite(var(check == 1))) > 0
var(check ~= 1) = Inf;
end
[~, i] = min(var);

idxClose = findClosePoints(x, y, K(i));
%}
idxSegm = ismember(1:length(x), idxClose);
end

function idxClose = findClosePoints(x, y, i)
dist = sqrt((x-x(i)*ones(size(x))).^2 + (y-y(i)*ones(size(y))).^2);
if i == 1
    minDist = dist(2);
elseif i == length(x)
    minDist = dist(end-1);
else
    minDist = min([dist(i+1), dist(i-1)]);
end
idxClose = find(dist < minDist);
i = 1;
while i < length(idxClose)
    idx = ismember(idxClose, idxClose(i)+1);
    idx = idx | ismember(idxClose, idxClose(i)-1);
    idxClose = idxClose(~idx);
    i = i + 1;
end



end

function [idxClose, iNew] = fixClose(idxClose, i, x, y, T)

if nargin < 5
T = round(median(diff(idxClose)));
else
    T = round(T);
end
if i > 0
iNew = idxClose(i) + T;
else 
iNew = idxClose(1) - T;
end
eps = round(T/100);

if iNew > 1 && iNew < length(x);
k1 = max(1, iNew-eps);
k2 = min(iNew + eps, length(x));
for k = k1:k2
distK(k - k1 +1) = sum(sqrt((x(idxClose)-x(k)*ones(size(idxClose))).^2 + (y(idxClose)-y(k)*ones(size(idxClose))).^2));
%dist(2) = sum(sqrt((x(idxClose)-x(i+T)*ones(size(idxClose))).^2 + (y(idxClose)-y(i+T)*ones(size(idxClose))).^2));
%dist(3) = sum(sqrt((x(idxClose)-x(i+T+1)*ones(size(idxClose))).^2 + (y(idxClose)-y(i+T+1)*ones(size(idxClose))).^2));
end

[~, k] = min(distK);

iNew = k + k1 - 1;
idxClose(end+1) = iNew;
idxClose = sort(idxClose);

else
    iNew = 0;
end

end

function [varr, check] = checkPoints(idxClose, m, T)

check = true;
segms = diff(idxClose);
varr = var(segms);
segmT = mean(segms);
nLoops  = m/segmT;

if length(idxClose) < 3
 varr = Inf;
end
%
if  sum(segms < T/2) > max(nLoops*0.1,1) || sqrt(varr) > segmT/2 || abs(length(idxClose)-nLoops) > max(nLoops*0.1,1)
 check = false;
end
%}
end


function plotPCs(PCs, nPCs, phi)

PCpairs = combntns(1:nPCs, 2);

if groupSize == 2
    for pair = PCpairs'
        if tan(phi) > 1
            x1 = min(PCs(pair(2),:))/tan(phi);
            x2 = max(PCs(pair(2),:))/tan(phi);
            y1 = min(PCs(pair(2),:)); 
            y2 = max(PCs(pair(2),:));
        else
            x1 = min(PCs(pair(2),:));
            x2 = max(PCs(pair(2),:));
            y1 = min(PCs(pair(2),:))*tan(phi); 
            y2 = max(PCs(pair(2),:))*tan(phi);
        end
        figure;
        plot(PCs(pair(1),:), PCs(pair(2),:), 'b-', 'linewidth', 2);
        hold on;
        plot([x1, x2], [y1, y2], 'r-', 'linewidth', 2);
        axis tight;
        xlabel('Principal component, 1', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('Principal component, 2', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
        set(gca, 'FontSize', 16, 'FontName', 'Times')
    end  
end
end


function [idxSegm, phi] = segment1(pc, phi)

if nargin < 2
    phi = findSymmetryPars(pc(1,:), pc(2,:));
elseif isempty(phi)
    phi = findSymmetryPars(pc(1,:), pc(2,:));
end
if phi > 1
    idxDown = pc(2,:)/tan(phi) < pc(1,:);
else
    idxDown = pc(2,:) < pc(1,:)*tan(phi);
end
idxEnd = (idxDown - [0 idxDown(1:end-1)]) == 1;
idxSegm = ~idxEnd;

end

function idxSegm = segment2(pc)
%

idx1 = sign(pc(2,:)./pc(1,:));
    idx2 = sign(pc(2,:));
    idx1 = diff(idx1) == 2;
    idx2 = diff(idx2) == 2;
    idxSegm = idx1 == idx2;
%}
     
end