
% do:
%       y = H*L*H*x
%       H is the Householder matrix
%
%   H = (I - u*u')
%   I = identity matrix
%   u = a vector, called the Householder vector
%   u*u' = an outer product
%
%   y = H*x
%   norm(y,2) == norm(x,2)
%   also y(1) is set to 0
%
%
%   H*L*H == L - u*v' - v*u'
%
%   H = (I - u*u'/alpha)
%   v = funky stuff

% in mypcg2:

q = myfunc (L,u,v,p)   

    q = (L*p) - u*(v'*p) - v*(u'*p) ;

% in hdip_fiedler: y = H*L*H*x

y = myfunc (L,u,v,x)

    y = (L*x) - u*(v'*x) - v*(u'*x) ;

% z = myfunc (H,C,r) to compute z = H*C*H*r

    y = happly (u, alpha, r) ;      % y = H*r
    y = C*y ;                       % y = C*y
    z = happly (u, alpha, y) ;      % z = H*y

%   H*L*H is a dense matrix, n^2 entries, n^2 >> # edges
%   L is sparse, with # entries = 2 * (# edges) << n^2
%   y = H*L*H*x
%   y = (L - u*v' - v*u')*x
%   y = (L*x) - u*(v'*x) - v*(u'*x) ;
%   e = nnz (L) = 2 * #edges
%
%   work:   L*x is O(e)
%           u*(v'*x) is O(n)
%           v*(u'*x) is O(n)
%           total is O(e+n)  <<< O(n^2) since e << n^2
%       

%   y = H*x = 
%           y = x - u * ((u'*x) / alpha) ;

%   y = H*L*H*x = H * (L * (H*x))
%
%       t = H*x, one call to happly:  2 dots/scale (like u*(u'*x))
%       s = L*t, one mxv
%       y = H*s, one call to happly:  2 dots/scale
%       total: 1 mxv, 4 dots/scale

%   y = (L*x) - u*(v'*x) - v*(u'*x) ;
%       1 mxv   dot         dot/scale
%       total: 1 mxv, 2 dots/scale




