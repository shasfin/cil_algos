function [U,z] = KMEANS(X,K)
% KMEANS(X,k) clusters the column of X into k clusters
% using the K-means algorithm
% The columns of U contain the centroids
% z is the assignment vector
tol = 1e-4;
change = 2*tol;
score_new = Inf;
max_it = 20;

[D,N] = size(X);

% initialize the centroids to random data points and preallocate the
% assignment vector
ind = randperm(N); % permute columns of X
U = X(:,ind(1:K)); % take the first k of the permuted columns
z = zeros(N,1);

it = 0;

while (it < max_it && change > tol)
    it = it + 1;
    % optimize Z (assign points to nearest cluster)
    Dist = repmat(sum(X.^2,1),K,1) - 2*U'*X + repmat(sum(U.^2,1)',1,N);
    [~,z] = min(Dist,[],1); 
    % optimize U (set centroids to the cluster means)
    for k=1:K
        U(:,k) = max(0,mean(X(:,z==k),2));
    end
    % update the score
    score_old = score_new;
    score_new = sum(sum((X - U(:,z)).^2));
    change = abs(score_old - score_new);
    disp(sprintf('Iterations = %d  Change = %0.5g.', it, change));
end