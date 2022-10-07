import pygraphblas as gb
import laplacian as lap
from norm2 import norm2
from math import sqrt
from happly import happly 
from hmhx import hmhx 

def mypcg2(L,u,malpha,invdiag,b,tol,maxit):
    '''
    mypcg2: Preconditioned conjugate gradient
    
    [x,k] = mypcg2 (L,u,v,alpha,C,b,tol,maxit)
    
    Solves A*x=b using pre-conditioned conjugage gradient and a diagonal
    preconditioner, where A = H*L*H and H is a Householder reflection with
    H = I - u*u'/alpha.  A = H*L*H is equal to L-u*v'-v*u' but this is not
    formed explicitly.  C is the pre-conditioner, and is a diagonal matrix with
    C(i,i) = 1/L(i,i).  Thus, H*C*H is an approximation of the "inverse" A.
    More precisely if F = H*C*H then F is an approximation inv (A (2:n,2:n)),
    because A itself is exactly singular.
    
    L is the Laplacian matrix of an undirected graph, with L(i,i) = the degree of
    node i, and L(i,j) = -1 for each edge (i,j).  L(i,i) must be > 0 for all i.
    
    Returns k as the # of iterations taken.
    
    From Golub and Van Loan, 2nd Edition, Algorithm 10.3.1, p529,
    
    Note that in the paper by Wu et al., the system A2*x=b has dimension n-1,
    with A2 = A (2:n,2:n).  Here, all of A is handled, but A (:,1) and A (1,:)
    are all zero, and so is b (1) and the solution x (1).  A2 must be symmetric
    and positive definite.
    
    In the GraphBLAS version, the input b can be overwritten with the solution x.
   
    Example:
    from operator import itemgetter
    J=sorted(list(gb.Matrix.ssget(2760)), key=itemgetter(0))[0]
    J[1].print(level=3)
    pattern=J[1].pattern(gb.types.FP32)
    h=pattern+pattern.transpose()
    omega=lap.laplacian(pattern)
    n=omega[0].nrows
    u=gb.Vector.dense(gb.types.FP32,n,fill=1)
    u[0]=1+sqrt(n)
    alpha=n+sqrt(n)
    c=omega[0].diag()
    indiag=(1/c.cast(gb.FP32)).diag()
    x=gb.Vector.dense(gb.types.FP32,n,fill=1)
    x[0]=0
    mypcg2=(omega[0],u,alpha,indiag,x,.000001,50)
    print(mypcg2)
    
    First we load in the matrix into python from suitesparse.
    Second we make sure that the matrix is of type FP32 (can be also be FP64)
    We make sure the matrix is symmetric, by adding it with its transpose
    Take the laplacian of that matrix.
    Create a vector filled with ones with the same size and number of rows in matrix
    Change first value in vector to (1+sqrt(num of rows in matrix))
    Create alpha, by (num of rows in matrix+sqrt(num of rows in matrix))
    Extract the diagonal of the modified matrix
    Create the inverse of that diagonal
    Create another vector filled with ones, with zero being its first value.
    Lastly call mypcg2 with the created variables as its parameters.. 

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
    '''
    print("---------------------------------------------------------------------------------------------------------------")
    print("ENTERING MYPCG2")
    print("---------------------------------------------------------------------------------------------------------------")
    #problem size
    n= L.nrows
    #b[0]=0 is required for the input
    b[0]=0
    #initial residual, assuming x=0 initially
    r=b
    #steper = initial guess of the solution (all zero)
    steper=gb.Vector.dense(gb.types.FP32,n,fill=0)
    #initial rho
    rho=1
    for k in range(1,maxit+1):
        print("---------------------------------------------------------------------------------------------------------------")
        print("ITERATIONS IN MYPCG2:"+str(k))
        print("---------------------------------------------------------------------------------------------------------------")
        #apply the preconditioner, using hmhx
        z=hmhx(invdiag,u,r,malpha)
        z[0]=0
        '''
        --------------------------------------------------------------------------
         determine the next search direction, p
        --------------------------------------------------------------------------
        '''
        #save the prior rho
        rho_prior=rho      
        #compute the new rho
        rho=r.emult(z).reduce_float(mon=gb.types.FP32.PLUS_MONOID) 
        if(k==1):
            #first step is the direction p=z
            p=z
        else:
            #subsequent step in direction p
            beta = rho/rho_prior
            p=z+beta*p
        '''
        --------------------------------------------------------------------------
         apply the matrix: q = A*p
        --------------------------------------------------------------------------
        L*p is a matrix-vector multiply and is the major work for each iteraion.
        v'*p and u'*p are both vector dot products, resulting in a scalar.
        '''
        #hmhx is used on q
        q=hmhx(L,u,p,malpha)
        q[0]=0
        '''
        --------------------------------------------------------------------------
         determine the stepsize
        --------------------------------------------------------------------------        
        '''
        #gamma =p'*A*p
        gamma=p.emult(q).reduce_float(mon=gb.types.FP32.PLUS_MONOID)
        #stepsize to take
        stepsize = rho/gamma
        '''
        --------------------------------------------------------------------------
        update x along the search direction p, and update the residual
        --------------------------------------------------------------------------
        '''
        #takes a step towards the solution
        steper=steper+stepsize*p
        #updates the residual
        r=r-stepsize*q
        #steper[0] and r[0] is always 0
        steper[0]=0
        r[0]=0
    
        '''
        --------------------------------------------------------------------------
        check for termination: when the residual norm(r) = norm(b-A*x) is small
        --------------------------------------------------------------------------
        '''
        rnorm=norm2(r)
        if (rnorm<tol):
            break
    print("---------------------------------------------------------------------------------------------------------------")
    print("ABOUT TO EXIT MYPCG2")
    print("---------------------------------------------------------------------------------------------------------------")
    return steper,k

