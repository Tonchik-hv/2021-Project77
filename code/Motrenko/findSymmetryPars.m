function alpha = findSymmetryPars(x, y)

[r, phi] = xy2rphi(x, y);
n = length(x);

R = [];
for alpha = 0:0.1:1
    xx = [];
    yy = [];
    minus = 0;
    [xAl, yAl] = rphi2xy(r, phi + alpha*pi/2);
    if  alpha ~= 0
        idxDwn = xAl*tan(pi/2 - alpha) > yAl;
        idxUp = xAl*tan(pi/2 - alpha) < yAl;
    else
        idxDwn = xAl > 0;
        idxUp = xAl < 0;
    end
    xx(1:sum(idxDwn)) = xAl(idxDwn);
    xx(end+1:end+sum(idxUp)) = xAl(idxUp);
    xx(end + 1:length(x)) = xAl(~idxDwn & ~idxUp);    
    yy(1:sum(idxDwn)) = yAl(idxDwn);
    yy(end+1:end+sum(idxUp)) = yAl(idxUp);
    yy(end + 1:length(x)) = yAl(~idxDwn & ~idxUp);
    
    m = sum(xAl ~= 0);
    if mod(m, 2) == 1
        minus = 1;
        xx = xx(2:end);
        yy = yy(2:end);
    end
    Q = constructQ(floor(m/2), n - minus);
    Xs = Q*[xx yy]';
    R(end+1) = sqrt(sum((Xs - [xx yy]').^2));
end

[~, idxMin] = min(R);
alpha = pi*(1 - 0.1*idxMin)/2;
end

function Q = constructQ(m,n)

S = zeros(m);
S(:, m+1:2*m) = eye(m);
S(m+1:2*m, 1:m) = eye(m);
S(2*m+1:n, 2*m+1:n) = eye(n - 2*m);


SS = -S;
SS(n+1:2*n, n+1:2*n) = S;
Q = 0.5*eye(2*n) + 0.5*SS;

end


function [r, phi] = xy2rphi(x, y)

r = sqrt(x.^2 + y.^2);
phi(x ~= 0) = atan(y(x~=0)./x(x~=0));
phi(x == 0) = pi/2;

end


function [x, y] = rphi2xy(r, phi)

x = r.*cos(phi);
y = r.*sin(phi);

end