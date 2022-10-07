function [x,lambda] = hd_orig (L)
%HDIP_FIEDLER: compute the Fiedler vector of a graph via the HDIP method.
% x = hdip_fiedler (L) computes the Fiedler vector of an undirected graph.
% The input matrix L is the Laplacian of the graph: a symmetric matrix
% with L(i,i) = the degree of node i in the graph (no self edges),
% and L(i,j) = L(j,i) = -1 if the edge (i,j) exists in the graph.
%
% The HDIP method is presented in the following paper:
%
%       ----------------------------------------------------------------------
%       Jian-ping Wu, Jun-qiang Song, Wei-min Zhang, An efficient and
%       accurate method to compute the Fiedler vector based on Householder
%       deflation and inverse power iteration, Journal of Computational and
%       Applied Mathematics, Volume 269, 2014, Pages 101-108, ISSN 0377-0427,
%       https://doi.org/10.1016/j.cam.2014.03.018.
%
%       Abstract: The Fiedler vector of a graph plays a vital role in many
%       applications. But it is usually very expensive in that it involves
%       the solution of an eigenvalue problem. In this paper, we introduce
%       the inverse power method incorporated with the Householder deflation
%       to compute the Fiedler Vector. In the inverse power iterations, the
%       coefficient matrix is formed implicitly, to take advantage of the
%       sparsity. The linear systems encountered at each iteration must be
%       symmetric positive definite, thus the conjugate gradient method is
%       used. In addition, preconditioning techniques are introduced to
%       reduce the computation cost. Any kind of preconditioning techniques
%       with dropping can be used. For the graphs related to some of the
%       sparse matrices downloaded from the UF Sparse Matrix Collection, the
%       experiments are compared to the known novel schemes. The results show
%       that the provided method is more accurate. While it is slower than
%       MC73 sequentially, it has good parallel efficiency compared with
%       TraceMin.  In addition, it is insensitive to the selection of
%       parameters, which is superior to the other two methods.
%       ----------------------------------------------------------------------
%
% This version is slow.  It is a direct translation of Figure 2.1 of the
% paper, but using brute force methods.  L2 should not be explictly formed,
% and L2\t is a direct method.  The paper uses preconditioned conjugate
% gradient.  See hdip_fielder for the method that has these features.
%
% The output vector x can also be obtained from the MATLAB eigs function:
%
%       L = laplacian (G) ;
%       [V,D] = eigs (L, 2, 'smallestabs') ;
%       x = V (:,2) ;
%
% Example:
%
%   A = bucky ;
%   G = graph (A) ;
%   figure (1) ; subplot (2,2,1) ; plot (G) ;
%   L = laplacian (G) ;
%   [x, lambda] = hd_orig (L) ;             % use HDIP to find Fiedler vector
%   [V,D] = eigs (L, 2, 'smallestabs') ;    % use eigs to find Fiedler vector
%   x2 = V (:,2) ;
%   lambda
%   lambda2 = D (2,2)
%   lambda - lambda2
%   % the sign of x and x2 can differ, so give x2 the same sign as x
%   x2 = x2 * sign (x (1)) ;
%   [x x2]
%   norm (L*x  - lambda*x, 2)
%   norm (L*x2 - lambda2*x2, 2)
%   % partition the graph
%   mid = median (x) ;
%   red = [1 0 0] ;
%   green = [0 1 0] ;
%   n = length (x) ;
%   color = zeros (n, 3) ;
%   left = find (x <= mid) ;
%   right = find (x > mid) ;
%   color (left, 1) = 1 ;
%   color (right, 2) = 1 ;
%   subplot (2,2,2) ; plot (G, 'NodeColor', color) ;
%   p = [left ; right] ;
%   S = A (p,p) ;
%   nleft = length (left) ;
%   edges_cut = nnz (S (1:nleft, nleft+1:n))
%   subplot (2,2,3) ; spy (A) ;
%   subplot (2,2,4) ; spy (S) ;

% Figure 2.1 of the paper gives the details of the method, translated below. 
% Line numbers 1. to 11. refer to the lines in Figure 2.1.

kmax = 30 ;
emax = 1e-6 ;

% 1. Set u(1) = 1+sqrt(n),  u(2:n) = 1, and alpha = n+sqrt(n)
n = size (L, 1) ;
u = ones (n, 1) ;
u (1) = 1 + sqrt (n) ;
alpha = n + sqrt (n) ;

% 2. Compute h = L*u/alpha, gamma = u'*h/alpha and v = h - gamma*u/2
h = L*(u/alpha) ;
gamma = (u'*h)/alpha ;
v = h - (gamma/2)*u ;

% 3. Set t(1:n-1) = 1
t = ones (n-1,1) ;

% compute the inf-norm of L:
Lnorm = norm (L, inf) ;

% this is just for illustration, L2 must not be explicitly computed:
r = u (2:n) ;
s = v (2:n) ;
L2 = L (2:n, 2:n) - r*s' - s*r' ;   

% 4. for k = 1:kmax
for k = 1:kmax

    % 5. compute beta = norm (t, 2) and set t = t/beta

    beta = norm (t, 2) ;
    t = t / beta ;

    % 6. compute y = L2*t, lambda = t'*y, and e = norm (y-lambda*t,inf)
    y = L2 * t ;
    lambda = t'*y ;
    e = norm (y - lambda*t, inf) ;

    % 7. if e / ||L||_inf < emax, break
    if (e / Lnorm < emax)
        break ;
    end

    % 8. solve L2*z = t for z, and set t = z
    z = L2 \ t ;    % this is a hack; need to use PCG here
    t = z ;

% 9. enddo
end

% 10. set x(1) = 0, "x(1:n) = t(1:n-1)" and compute beta = u'*x/alpha
% (typo in the paper: should be "x(2:n) = t(1:n-1)")
x = [0 ; t] ;
beta = (u'*x) / alpha ;

% 11. set x = x - beta*u
x = x - beta*u ;

