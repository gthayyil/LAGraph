% to use type these MATLAB commands:

%   vi ~/Documents/MATLAB/startup.m
%   with this line:
%
%       addpath ~/ssget

%   diary mystuff.txt
%   domats

clear all

index = ssget ;
list = find (index.pattern_symmetry == 1 & index.nnz > 1e6) ;
% & index.ncc == 1) ;
[ignore, i] = sort (index.nnz (list)) ;
list = list (i)' ;

list = [1883 1361]
nmatrices = length (list) ;

for k = 1:nmatrices
    id = list (k) ;
    fprintf ('%4d %d %s/%s\n', id, index.nnz (id), index.Group{id}, ...
        index.Name{id}) ;
    Prob = ssget (id)
    A = Prob.A ;
    A = spones (A) ;
    if (isfield (Prob, 'Zeros')) 
        Z = Prob.Zeros ;
        A = A + Z ;
    end
    A = tril (A, -1) ;
    A = A+A' ;
    fprintf ('nvals(A): %d\n', nnz (A)) ;

    G = graph (A) ;
    L = laplacian (G) ;
    % spy (L)
    t1 = tic ;
    % use HDIP to find Fiedler vector
    [x, lambda, iters] = hdip_fiedler (L, [ ], 1e-6, 1e-6) ;
    t_hdip_fiedler = toc (t1) ;
    fprintf ('lambda: %g iters: %d %d hdip_fiedler time: %g\n', lambda, iters, t_hdip_fielder) ;

    t1 = tic ;
    [ignore, p] = sort (x) ;
    S = A (p,p) ;
    n = size (S,1) ;
    nleft = floor (n/2) ;
    left = 1:nleft ;
    right = (nleft+1):n ;
    edge_cut = nnz (S (left, right)) ;
    tcut = toc (t1)

    fprintf ('edge cut: %d cut time: %g\n', edge_cut, tcut) ;

end
