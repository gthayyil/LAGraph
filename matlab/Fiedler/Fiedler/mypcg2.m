function [x,k] = mypcg2 (L,u,v,alpha,C,b,tol,maxit)
%MYPCG2:  preconditioned conjugate gradient
%
%   [x,k] = mypcg2 (L,u,v,alpha,C,b,tol,maxit)
%
% Solves A*x=b using pre-conditioned conjugage gradient and a diagonal
% preconditioner, where A = H*L*H and H is a Householder reflection with
% H = I - u*u'/alpha.  A = H*L*H is equal to L-u*v'-v*u' but this is not
% formed explicitly.  C is the pre-conditioner, and is a diagonal matrix with
% C(i,i) = 1/L(i,i).  Thus, H*C*H is an approximation of the "inverse" A.
% More precisely if F = H*C*H then F is an approximation inv (A (2:n,2:n)),
% because A itself is exactly singular.
%
% L is the Laplacian matrix of an undirected graph, with L(i,i) = the degree of
% node i, and L(i,j) = -1 for each edge (i,j).  L(i,i) must be > 0 for all i.
%
% Returns k as the # of iterations taken.
%
% From Golub and Van Loan, 2nd Edition, Algorithm 10.3.1, p529,
%
% Note that in the paper by Wu et al., the system A2*x=b has dimension n-1,
% with A2 = A (2:n,2:n).  Here, all of A is handled, but A (:,1) and A (1,:)
% are all zero, and so is b (1) and the solution x (1).  A2 must be symmetric
% and positive definite.
%
% In the GraphBLAS version, the input b can be overwritten with the solution x.

n = size (L,1) ;                    % problem size

% r is the residual, r = b-A*x, for the current solution x:
assert (b (1) == 0) ;               % requre b(1)=0 on input
r = b ;                             % initial residual, assuming x = 0 initially

% x = initial guess of the solution (all zero):
% note that at this point b is not needed, and x can overwrite it
x = zeros (n,1) ;                   % initial guess of solution, x = 0

rho = 1 ;                           % initial rho

for k = 1:maxit

    %--------------------------------------------------------------------------
    % apply the preconditioner: z = H*C*H*r
    %--------------------------------------------------------------------------

    % 1 mxv, 4 dot/scale ops
    % z = myfunc (H,C,r) to compute z = H*C*H*r
    y = happly (u, alpha, r) ;      % y = H*r
    y = C*y ;                       % y = C*y
    z = happly (u, alpha, y) ;      % z = H*y
    z (1) = 0 ;                     % z(1) may have roundoff noise, but is
                                    % algebraicly exactly zero

    %--------------------------------------------------------------------------
    % determine the next search direction, p
    %--------------------------------------------------------------------------

    rho_prior = rho ;               % save the prior rho
    rho = r'*z ;                    % compute the new rho
    if (k == 1)
        p = z ;                     % first step is the direction p = z
    else
        beta = rho / rho_prior ;    % subsequent step in direction p
        p = z + beta * p ;
    end
    assert (p (1) == 0) ;

    %--------------------------------------------------------------------------
    % apply the matrix: q = A*p
    %--------------------------------------------------------------------------

    % L*p is a matrix-vector multiply and is the major work for each iteraion.
    % v'*p and u'*p are both vector dot products, resulting in a scalar.

    % q = myfunc (H,L,p) to compute q = H*L*H*p
    q = (L*p) - u*(v'*p) - v*(u'*p) ;
    q (1) = 0 ;

    %--------------------------------------------------------------------------
    % determine the stepsize
    %--------------------------------------------------------------------------

    gamma = p' * q ;                % gamma = p'*A*p
    stepsize = rho / gamma ;        % stepsize to take 

    %--------------------------------------------------------------------------
    % update x along the search direction p, and update the residual
    %--------------------------------------------------------------------------

    x = x + stepsize * p ;          % take a step towards the solution
    r = r - stepsize * q ;          % update the residual
    assert (x (1) == 0) ;           % x(1) and r(1) always stay zero
    assert (r (1) == 0) ;

    %--------------------------------------------------------------------------
    % check for termination: when the residual norm(r) = norm(b-A*x) is small
    %--------------------------------------------------------------------------

    rnorm = norm (r, 2) ;
    if (rnorm < tol)
        break ;
    end
end

