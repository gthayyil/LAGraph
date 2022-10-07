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
    tic
    % use HDIP to find Fiedler vector
    [x, lambda, iters] = hdip_fiedler (L, [ ], 1e-6, 1e-6) ;
    toc
    fprintf ('lambda: %g iters: %d %d\n', lambda, iters) ;

    [ignore, p] = sort (x) ;
    S = A (p,p) ;
    n = size (S,1) ;
    nleft = floor (n/2) ;
    left = 1:nleft ;
    right = (nleft+1):n ;
    edge_cut = nnz (S (left, right)) ;
    fprintf ('edge cut: %d\n', edge_cut) ;

end
