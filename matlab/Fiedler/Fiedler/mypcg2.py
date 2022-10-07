
import pygraphblas as gb
import laplacian as lap
from norm2 import norm2
from math import sqrt
from happly import happly
from hmhx import hmhx
def mypcg2(L,u,malpha,invdiag,b,tol=.000001,maxit=50):
    n= L.nrows
    #print(n)
    b[0]=0
    r=b
    steper=gb.Vector.dense(gb.types.FP64,n,fill=0)
    rho=1
    for k in range(1,maxit+1):
        z=hmhx(invdiag,u,r,malpha)
        z[0]=0
        #test=happly(u,r,malpha)
        #print(test)
        #test=invdiag.mxv(test)
        #print(test)
        #print(invdiag)
        #print(z)
        #wait = input("Press Enter to continue.")

        rho_prior=rho
        rho=r.emult(z).reduce_float()
        if(k==1):
            p=z
        else:
            #print(p)
            #print(z)
            beta = rho/rho_prior
            p=z+beta*p

        q=hmhx(L,u,p,malpha)
        q[0]=0
        #print(q)
        #print(p)
        #wait = input("Press Enter to continue.")
        gamma=p.emult(q).reduce_float()
        #print(gamma)
        #print(rho)
        #wait = input("Press Enter to continue.")
        stepsize = rho/gamma
        #print(stepsize)
        #print(steper)
        #print()
        #print(stepsize)
        steper=steper+stepsize*p
        r=r-stepsize*q
        steper[0]=0
        r[0]=0
        #print(steper)
        #print(r)
        #wait = input("Press Enter to continue.")

        rnorm=norm2(r)
        if (rnorm<tol):
            break
    return steper,k

