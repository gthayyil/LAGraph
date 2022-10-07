index = ssget ;
list = find (index.pattern_symmetry == 1 & index.ncc == 1 & index.nnz>1e6 ) ;
[ignore, i] = sort (index.nnz (list)) ;
list = list (i)' ;
nmatrices = length (list)

list = [2516 1883 1882 1899 1361 1893 1861 1454 1385 1257]
nmatrices = length (list) ;

for k = 1:nmatrices %143:nmatrices
	id = list (k) ;
	id
	Prob = ssget (id) ;
	A = Prob.A ;

	% A = Prob.A ;
	A=spones(A);
	A=A+A';
	A=spones(A);
	A=tril(A,-1);
	A=A+A';

	t0 = tic ;
	mystuff(A);
	t = toc (t0) ;
	fprintf('The following time is the total time: %g\n', t)
	 pause;
end
