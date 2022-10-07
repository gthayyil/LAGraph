function y = happly (u, alpha, x)
%HAPPLY apply a Householder reflection
% y = happly (u, alpha, x)
% computes y = H*x where H = (I - u*u'/alpha)

% H*x = (I - u*u'/alpha) * x = (x - u*(u'*x) / alpha)

% y, u, and x are vectors of size n, and alpha is a scalar
% in the GraphBLAS version, the inputs u, alpha, and x must not
% be modified, so y must be a different vector.

y = x - u * ((u'*x) / alpha) ;

%   s = sum (x) is almost u'*x
%   s = ((u'*x) / alpha) ;
%
%   assign scalar to x with accum:
%   x -= scalar 
%
%       does:
%
%           for i = 1:n
%               x (i) = x(i) - scalar


