function X = hankelmatrix(x,k)
% h = hankelmatrix(H,K) converts column-vector to Hankel matrix and back again, 
% SSA-modification averages matrix anti-diagonals to make a vector
%
% x [m,1] vector (time series of SSA)
% X [L,K] where L = N-K+1 hankel matrix
% K [int scalar] number of cloumns in the matrix
%
% Example
% x = [1:17];
% k = 5;
% X = hankelmatrix(x,k);
% x1 = hankelmatrix(X);
%
% http://strijov.com 05-Jul-2008

%1) from vector to matrix
[L,K] = size(x);
if min(L,K) == 1 % 1) from vector to matrix
    N = length(x);
    L = N-k+1;
    X = NaN*ones(L,k); % returns NaNs if h contains NaNs
    for i=1:L % WARNING the history begins at 1-st sample of the ts
        X(i,:) = x(i:i+(k-1));
    end
else     
% 2) from matrix to vector and average anti-diagonals
    N = L+K-1;
    X = zeros(N,K); % now X will be time series and x is the Hankel matrix
    % input argument k is not used
    for j=1:K
        X(j:j+L-1,j) = x(:,j);
    end
    denom = [1:K K*ones(1,N-2*K) K:-1:1]'; % number of non-zero (non-empty) items in the matrix
    % TODO denom works only for L>K 
    X = sum(X,2)./denom; % SSA diagonal averaging
end%if

