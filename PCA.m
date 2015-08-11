function Xapprox = PCA(X,k)
figure
imshow(X)

% approximate data matrix X using only k dimensions with PCA
[~,n] = size(X);

% build Xbar
Xbar = X - repmat(mean(X),n,1);

% build the covariance matrix
CovX = Xbar * Xbar' / n;

% retrieve the eigendecomposition
[U,D] = eig(CovX);

% sort in descending order
[~,I] = sort(diag(D),'descend');
U = U(:,I);

% approximate X
Xapprox = U(:,1:k)*(U(:,1:k)' * Xbar) + repmat(mean(X),n,1);

figure
imshow(Xapprox)

end