
% To compute:

    beta = norm (x, 2) ;

% do this in GraphBLAS:

    beta = sqrt (sum (x.^2)) ;

% which is done like this:

    q = x.^2 ;  % q.apply_2nd with POW op, scalar 2
    % q(i) is now pow(x(i),2) for all i

    ss = sum(q) % q.reduce with PLUS monoid to get scalar ss
                % (use MAX for inf-norm)

    % sqrt of the scalar:
    beta = sqrt (ss) ;

% to compute the inf-norm:

    e = norm (y - lambda*x, inf) ;

% do this:

    q = abs (y - lambda*x)
    e = max (q)

