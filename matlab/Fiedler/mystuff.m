function[iters,lambda,timerVal,S,edges_cut]=mystuff(A) 
 % Prob = ssget(id);
 tstart1 = tic ;
 tstart = tic ;

 G=graph(A);
 L = laplacian (G);
 tlap = toc (tstart) ;
 fprintf ('compute Laplacian: %g sec\n', tlap) ;
 tstart = tic ;
 [x, lambda, iters] = hdip_fiedler (L);
 timerVal=toc (tstart);
 tstart = tic ;
 %{
 [V,D] = eigs (L, 2, 'smallestabs') ;    % use eigs to find Fiedler vector
 x2 = V (:,2) ;
 lambda2 = D (2,2);
 lambda - lambda2;
 % the sign of x and x2 can differ, so give x2 the same sign as x
 x2 = x2 * sign (x (1)) ;
 % [x x2]
 norm (L*x  - lambda*x, 2);
 norm (L*x2 - lambda2*x2, 2);
 % partition the graph
 %mid = median (x) 
 red = [1 0 0] ;
 green = [0 1 0] ;
 n = length (x) ;
 color = zeros (n, 3) ;
 %left = find (x <= mid) ;
 %right = find (x > mid) ;
 %}
 [y,p]=sort(x);
tsort=toc(tstart)
 n = length (x) ;
 half = floor (n/2) ;
 left = p(1:half) ;
 right = p((half+1):n) ;
 % color (left, 1) = 1 ;
 % color (right, 2) = 1 ;
 %% subplot (2,2,2) ; plot (G, 'NodeColor', color) ;
 % p = [left ; right] 
 S = A (p,p); 
 nleft = length (left); 
 edges_cut = nnz (S (1:nleft, nleft+1:n))
 tcut=toc (tstart);
 fprintf ('cut time: %g\n', tcut) ;
 fprintf ('sort time: %g\n', tsort) ;
 ttot=toc (tstart1);
 fprintf('The follwing time is just for hdip function: ')
 timerVal
 iters
 lambda
 fprintf('The following time is calculated without ssget: %g\n', ttot)
 % subplot (2,2,3) ; spy (A) ;
 % subplot (2,2,4) ; spy (S) ;
end
	      

