A = bucky ;
% A = A (1:16,1:16) ;
G = graph (A) ;
%% figure (1) ; subplot (2,2,1) ; plot (G) ;
L = laplacian (G) ;
[x, lambda, iters] = hdip_fiedler (L) ; % use HDIP to find Fiedler vector
iters
[V,D] = eigs (L, 2, 'smallestabs') ;    % use eigs to find Fiedler vector
x2 = V (:,2) ;
lambda
lambda2 = D (2,2)
lambda - lambda2
% the sign of x and x2 can differ, so give x2 the same sign as x
x2 = x2 * sign (x (1)) ;
% [x x2]
norm (L*x  - lambda*x, 2)
norm (L*x2 - lambda2*x2, 2)
% partition the graph
mid = median (x) ;
red = [1 0 0] ;
green = [0 1 0] ;
n = length (x) ;
color = zeros (n, 3) ;
left = find (x <= mid) ;
right = find (x > mid) ;
color (left, 1) = 1 ;
color (right, 2) = 1 ;
%% subplot (2,2,2) ; plot (G, 'NodeColor', color) ;
p = [left ; right] ;
S = A (p,p) ;
nleft = length (left) ;
edges_cut = nnz (S (1:nleft, nleft+1:n))
%% subplot (2,2,3) ; spy (A) ;
%% subplot (2,2,4) ; spy (S) ;

