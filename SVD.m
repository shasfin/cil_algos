function X_approx = SVD(X)
%use svd decomposition to build a low-rank approximation of X
% infer k from the fraction to the first singular value

r = 0.1; % keep only those singular values which are at most r times bigger than the largest singular value

[U,D,V] = svd(X,'econ');
d = diag(D);
d = d / d(1);
k = sum(d > r); % infer k

X_approx = U(:,1:k)*D(1:k,1:k)*V(:,1:k)';

end

