
import pygraphblas as gb
#from operator import itemgetter
def laplacian(G):
    A=G.offdiag()
    d=gb.Vector.dense(gb.types.FP64,A.ncols)
    #print(d)
    A.mxv(d,accum=gb.types.FP64.PLUS,semiring=gb.types.FP64.PLUS_PAIR,out=d)
    #print("This is d after A.mxv(d..)")
    #print(d)
    S=gb.Matrix.sparse(gb.types.FP64,A.nrows,A.ncols)
    S.assign_scalar(-1,mask=A,desc=gb.descriptor.S)

    # d [mask is d, complemented: d==0] = 1

    Dmatrix=gb.Matrix.sparse(gb.types.FP64,A.nrows,A.ncols)
    inform=d.reduce_float(d.type.MAX_MONOID)*2
    #print("The inf norm is "+str(inform))
    for i in range(A.nrows):
        Dmatrix[i,i]=d[i]

#   I = [0, 1, ... n-1]
#   Dmatrix = build (0:n-1, 0:n-1, d)

#   v5.1
#   Dmatrix = diag (d,0)

    L=Dmatrix+S
    #print(A)
    #print(d)
    #print(S)
    #print(L)
    return L
def inform(A):
    A=A.offdiag()
    d=gb.Vector.dense(gb.types.FP64,A.ncols)
    #print(d)
    A.mxv(d,accum=gb.types.FP64.PLUS,semiring=gb.types.FP64.PLUS_PAIR,out=d)
    inform=d.reduce_float(d.type.MAX_MONOID)*2
    #print("The inf norm is "+str(inform))
    return inform
#J=sorted(list(gb.Matrix.ssget(2760)), key=itemgetter(0))[0]
#J[1].print(level=3)
#L=laplacian(J[1])
