function [L,S] = RPCA(X)
%ADMM method for RPCA
figure
imshow(X);

max_it = 100;
delta = 1e-7;

[m,n] = size(X);
rho = 5;  % from the paper: rho = (m*n/4) / norm(X,1);
lambda = 1 / sqrt(n);

S = zeros(size(X));
L = zeros(size(X));
nu = zeros(size(X));
t = 0;


while t <= max_it % from the paper norm(X-L-S,'fro') > delta * norm(X,'fro')
   
   L = shrink_sv(X - S - 1/rho * nu, 1/rho);
   S = shrink(X - L - 1/rho * nu, lambda / rho);
   nu = nu + rho*(L + S - X);
   
   t = t+1;
   
end

figure
imshow(L);
figure
imshow(S);
figure
imshow(L+S);


end

function Y = shrink(X,threshold)

    Y = sign(X) .* max(abs(X)-threshold,0);
    
end

function Y = shrink_sv(X,threshold)

[U,D,V] = svd(X);
Y = U * shrink(D,threshold) * V';

end

