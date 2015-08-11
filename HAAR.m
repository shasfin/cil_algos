function X_rec = HAAR(X)
%apply Haar transform (computed by another function) to an input matrix X
%threshold the coefficients
%transform back to normal domain

% assume X is quadratic, one channel, dimension is a power of two

figure
imshow(X);

r = 0.1; % relative threshold (keep a portion r of the values)
[m,n] = size(X);
H = haarTrans(m);

% apply H to rows and then to the columns of the resulting matrix
T = H' * X * H;

% threshold the coefficients
[~,ind] = sort(abs(T(:)),'descend');
T(ind(floor(r*length(ind)):end)) = 0;

% reconstruct image
X_rec = H*T*H';

figure
imshow(X_rec);

end