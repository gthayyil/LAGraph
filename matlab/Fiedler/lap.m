function L = lap(A)

d[i]


L = Matrix.sparse(types.FP64, n, n)
for i in [0:n-1]
    L[i,i] = d[i]

L.wait()

print(L)
