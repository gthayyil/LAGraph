
A = sparse (13, 13) ;
A (1:5, 1:5) = rand (5) ;
A (6:10, 6:10) = rand (5) ;
A (11:13, 11:13) = rand (3) ;
A = A+A' ;
A = tril (A,-1) + triu (A, 1) ;
A (:,1) = 0 ;
A (1,:) = 0 ;
spy (A)
G = graph (A) ;
L = laplacian (G)
[x, lambda, iters] = hdip_fiedler (L) ; % use HDIP to find Fiedler vector
iters
x
lambda
norm (L*x - lambda*x, 2)

[V,D] = eigs (L, 2, 'smallestabs') ;    % use eigs to find Fiedler vector
x2 = V (:,2)
lambda2 = D (2,2)
norm (L*x2 - lambda2*x2, 2)
