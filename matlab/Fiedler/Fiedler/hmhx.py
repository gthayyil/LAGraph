
from happly import happly
from norm2 import norm2
import pygraphblas as gb
def hmhx(C,u,r,alpha):
    y=happly(u,r,alpha)
    g=C.mxv(y)
    z=happly(u,g,alpha)
    z[0]=0
    #print(z)
    return z
I=[4.0,5.0,6.0,7.0]
v=gb.Vector.from_list(I)
x=gb.Vector.from_list(I)
#print(v)
A = 2.0
I=[3,3,3]
J=[1,2,3]
K=[4,5,6]
#x=gb.Vector.dense(gb.types.FP64,3,fill=1)
#x[0]=0
#print(x)
#beta=norm2(x)
#x=x/beta
#print(x)
M=gb.Matrix.from_lists(I,J,K)
M=M+M.transpose()
mdiag=M.diag()
indiag=1/mdiag.cast(gb.FP64)
#print(indiag)
#print(hmhx(indiag,v,x,A))