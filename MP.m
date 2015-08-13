function Z = MP(X)
%MP Matching Pursuit
%   Detailed explanation goes here
K = 5; % for the basket image, 64 is the largest value computed within few seconds
tol = 1e-4;

L = 5; % overcompletness factor
[D,N] = size(X);

U = overDCTdict(D,L*D);

R = X;
Z = zeros(L*D,N);

% iterate: find vector in U with largest scalar product with respect to R,
% adjust Residuum R and Code Z
while (sum(sum(Z~=0) < K) > 0 && norm(X-U*Z) >= tol )% stop when the zero-norm of each column of Z is >= K
    % find best atoms
    SP = U' * R;
    [maxvec,ind] = max(abs(SP));
    
    % update Z
    vals = sub2ind(size(Z), ind, 1:N);
    Z(vals) = Z(vals) + SP(vals);
    
    % update R
    R = R - repmat(maxvec,D,1).*U(:,ind);
    
end

end

