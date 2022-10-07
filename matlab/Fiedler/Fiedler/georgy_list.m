% create a list of symmetric matrices with one component in them

clear all
diary georgy_list.txt
type georgy_list.m
index = ssget
list = find (index.pattern_symmetry == 1 & index.ncc == 1) ;
[ignore, i] = sort (index.nnz (list)) ;
list = list (i)'
nmatrices = length (list)
for k = 1:nmatrices
    id = list (k) ;
    fprintf ('%4d %d %s/%s\n', id, index.nnz (id), index.Group{id}, ...
        index.Name{id}) ;
end


diary off
