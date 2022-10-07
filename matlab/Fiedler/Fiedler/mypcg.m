function [x,k] = mypcg (A,b,tol,maxit)
%MYPCG:  preconditioned conjugate gradient
%
%   [x,k] = mypcg (A,b,tol,maxit)
%
% Solves A*x=b using pre-conditioned conjugage gradient
% and a diagonal preconditioner.  A must be symmetric
% and positive definite.  tol defaults to 1e-6 if not present
% or if empty.  maxit defaults to 50 if not present or empty.
% Returns k as the # of iterations taken.
%
% From Golub and Van Loan, 2nd Edition, Algorithm 10.3.1, p529.
%
% Example:
%
%   n = 10 ;
%   A = rand (n) ;
%   A = A'*A ;
%   b = rand (n,1) ;
%   x1 = A\b ;
%   x2 = mypcg (A,b) ;
%   err = norm (x1-x2)
%   resid = norm (A*x2-b)


% tol = stopping criterion (default 1e-6)
if (nargin < 3)
    tol = [ ] ;
end
if (isempty (tol))
    tol = 1e-6 ;
end

% maxit = max number of iterations (default 30)
if (nargin < 4)
    maxit = [ ] ;
end
if (isempty (maxit))
    maxit = 50 ;
end

n = size (A,1) ;

% C = diagonal matrix with C(i,i) = 1/A(i,i)
c = diag (A) ;
if (any (c <= 0))
    error ('diag(A) must be positive') ;
end
C = diag (1./c) ;

x = zeros (n,1) ;                   % initial guess of solution
r = b ;                             % initial residual
rho = 1 ;                           % initial rho

for k = 1:maxit
    rho_prior = rho ;               % save the prior rho
    z = C * r ;                     % apply the preconditioner
    rho = r'*z ;                    % compute the new rho
    if (k == 1)
        p = z ;                     % first step in direction p = z
    else
        beta = rho / rho_prior ;    % subsequent step in direction p
        p = z + beta * p ;
    end
    q = A*p ;                       % apply the matrix
    ptAp = p' * q ;                 % ptAp = p'*A*p
    alpha = rho / ptAp ;            % alpha = stepsize to take 
    x = x + alpha * p ;             % take a step
    r = r - alpha * q ;             % update the residual
    rnorm = norm (A*x-b,2) ;        % check for termination
    if (rnorm < tol)
        break ;
    end
end

