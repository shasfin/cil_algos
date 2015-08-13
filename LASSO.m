function [x1,x2] = LASSO(A,b)
%ADMM for LASSO
[m,n] = size(A);

max_it = 200;
lambda = 0.5;
rho = 0.5;

x1 = zeros(n,1);
x2 = zeros(n,1);
nu = zeros(n,1);

At = A';
Atb = At*b;
[L,R] = lu(At*A + rho*eye(size(At*A)));


t = 0;
while t < max_it
    z = L \ (Atb + rho * x2 - nu);
    x1 = R \ z;
    
    x2 = shrink(x1 + nu / rho,lambda / rho);
    
    nu = nu + rho * (x1 - x2);
    
    t = t+1;
end

end

function Y = shrink(X,threshold)

    Y = sign(X) .* max(abs(X) - threshold,0);
    
end

