function x = ssaMtick(x,M,K,r)
% x = ssaMtick(x,M,K,r) run 1-tick SSA forecast function 
% 
% M [int scalar] number on the ticks to forecast after the end of the time series x
% x, K, r  see notations and example in ssa1tick
% 
if nargin < 3, error('not enough arguments'); end
if nargin == 3
    for m = 1:M;
        x = ssa1tick(x,K);
    end
else
    for m = 1:M;
        x = ssa1tick(x,K,r);
    end
end
    
return


function [x1, Lambda] = ssa1tick(x,K,r)
% [x1, Lambda] = ssa1tick(x,K,r) one-tick SSA forecast
% after Golyandina N. Metod Gusenitsa-SSA: Prognoz vremennych ryadov. S-Pb. 2004. P. 23.
% 
% x [N,1] time series, one-dimensional
% K [int scalar] period
% r [int scalar] rank of Hankel matrix 
%
% Example
% x = [1 2 3 2 1 2 3 2 1 2 3]';
% x = [0 1 2 2 1 0 0 1 2 2 1 0 0 1 2]';
% x1 = ssa1tick(x,4);
% t = 1:length(x)+1;
% plot(t',[x1 [x;NaN]]); hold on
% plot(t',[x1 [x;NaN]], '.'); hold off
% axis tight

if nargin < 2, error('more arguments required'); end
N = length(x);              % 
L = N-K+1;

X = hankelmatrix(x,K);      % make antidiagonal matrix
                            % U-columns are eigenvectors of XX', since
[U,Lambda,V] = svd(X);      % XX' = ULV' * VLU' = UL2U => XX'U = UL2
%idxr = find(diag(Lambda)>=Wmin);
idxr = find(diag(Lambda)>0);     % diag elements always decrease
r1 = max(idxr);
if nargin == 2 || r>r1, r=r1; end % align the number of the eigenvectors

X1 = U(:,1:r) * Lambda(1:r, 1:r) * V(:,1:r)'; % reconstruct the hankel matrix
x1 = hankelmatrix(X1);      % convert to the time series

pi = U(end,1:r);            % get the last component of the eigenvectors
Up = U(1:end-1,1:r);        % get the firt L-1 components of the eigenvectors

R = 1/(1-sumsqr(pi)) * (Up*pi'); % compute the recursive ts coeffiftients
%g = R'*x1(end-L+2:end);     % forecast os the weighted sum of the reconstructed ts
g = R'*x(end-L+2:end);     % try this (the raw ts)

x1 = [x1;g];                % append the forecast to the reconstructed ts
return
