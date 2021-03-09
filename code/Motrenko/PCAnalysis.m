function [PCs, ts_reduced, lstLambda] = PCAnalysis(ts, T)

[PCs, matLambda, ts_reduced] = PCA(ts', T, 1, 1);
lstLambda = diag(matLambda);

end

%
function [pc, Lambda, ts] = PCA(ts,K, ifCenter, ifNormal)
% [x1, Lambda] = ssa1tick(x,K,r) one-tick SSA forecast
% after Golyandina N. Metod Gusenitsa-SSA: Prognoz vremennych ryadov. S-Pb. 2004. P. 23.
% 
% x [N,1] time series, one-dimensional
% K [int scalar] period
% r [int scalar] rank of Hankel matrix 
%
% Example
% x = sin([1:48]*pi/6);
% x = [1 2 3 2 1 2 3 2 1 2 3]';
% x = [0 1 2 2 1 0 0 1 2 2 1 0 0 1 2]';
% x1 = ssa1tick(x,4);
% t = 1:length(x)+1;
% plot(t',[x1 [x;NaN]]); hold on
% plot(t',[x1 [x;NaN]], '.'); hold off
% axis tight

if nargin < 2, error('more arguments required'); 
end
if nargin < 4
    ifNormal = 0;
    ifCenter = 0;

end
N = length(ts);              % 
L = N-K+1;

X = [];
for d = 1:size(ts,1)
X = [X hankelmatrix(ts(d,:),K)];      % make antidiagonal matrix
end                            % U-columns are eigenvectors of XX', since
 
mean = [];
if ifCenter
    mean = sum(X', 2)/L;
    X = X-repmat(mean',L, 1);
end
err = [];
if ifNormal
    err = sqrt(sum(X'.^2, 2)/L);
    X = X./repmat(err', L, 1);
end                            
                            
[U,Lambda,V] = svd(X'*X/L, 0);      % XX' = ULV' * VLU' = UL2U => XX'U = UL2
%idxr = find(diag(Lambda)>=Wmin);
idxr = find(diag(Lambda)>10^(-6));     % diag elements always decrease
[Er, r] = findR(Lambda, 0.95);
Lambda(r+1:end, r+1:end) = zeros(length(Lambda)-r);
%r1 = max(idxr);
%if nargin == 2 || r>r1, r=r1; end % align the number of the eigenvectors


pc = U'*X';
pc = pc(idxr,:);

X = zeros(size(X));

for i = 1:r
    X = X + pc(i,:)'*U(:,i)';
end
[L, K] = size(X);
if L <= K
    disp('PCAanalysis: OOps!')
end
ts = hankelmatrix(X);

end



function plotLambda(lstLambda, nLam)

figure;
plot(sqrt(lstLambda(1:nLam)), 'b-', 'linewidth', 1.5);
hold on;
plot(sqrt(lstLambda(1:nLam)), 'k.', 'markersize', 10);
hold off;
axis tight;
xlabel('Principal component number, $j$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
ylabel('Sq. root of eigenvalue, $\sqrt{\lambda_j}$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter','latex');
set(gca, 'FontSize', 16, 'FontName', 'Times')

end

function [E r] = findR(Lambda, tau)

for r = 1:size(Lambda,2)
    E(r) = sum(diag(Lambda(1:r,1:r).^2))/sum((diag(Lambda)).^2);
end

idx = E >= tau;
r = min(find((1:size(Lambda,2)).*idx));

end

%{
function [tsMax, idxMax] = findMax(ts)

idxUp = (ts - [ts(1) ts(1:end-1)]) > 0;
idxDown = (ts - [ts(2:end) ts(end)]) > 0;
idxMax = idxUp & idxDown;
tsMax = ts(idxMax);

end
%}