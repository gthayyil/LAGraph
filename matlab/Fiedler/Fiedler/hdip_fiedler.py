import pygraphblas as gb
import laplacian as lap
from norm2 import norm2
from math import sqrt
from hmhx import hmhx
from mypcg2 import mypcg2
from operator import itemgetter
def hdip_fiedler(L,kmax =[20,50],emax = .000001,tol= .000001 ):
    """
    Hdip_fiedler: This function computes the Fiedler Vector of a graph by using the HDIP methos.
    x= hdip_fiedler(L) computes the Fiedler vector of the laplacian of a symmetric matrix.

    However, the laplacian function used here (created by the same authors) returns the
    laplacian of the matrix and its inf norm, in a tuple, as such the first input parameter
    to hdip_fiedler is a tuple with the laplacian of a symmetric matrix and its inf norm.


    The HDIP method used, was presented in this paper:
     ----------------------------------------------------------------------
            Jian-ping Wu, Jun-qiang Song, Wei-min Zhang, An efficient and
            accurate method to compute the Fiedler vector based on Householder
            deflation and inverse power iteration, Journal of Computational and
            Applied Mathematics, Volume 269, 2014, Pages 101-108, ISSN 0377-0427,
            https://doi.org/10.1016/j.cam.2014.03.018.

            Abstract: The Fiedler vector of a graph plays a vital role in many
            applications. But it is usually very expensive in that it involves
            the solution of an eigenvalue problem. In this paper, we introduce
            the inverse power method incorporated with the Householder deflation
 	    to compute the Fiedler Vector. In the inverse power iterations, the
            coefficient matrix is formed implicitly, to take advantage of the
            sparsity. The linear systems encountered at each iteration must be
            symmetric positive definite, thus the conjugate gradient method is
            used. In addition, preconditioning techniques are introduced to
            reduce the computation cost. Any kind of preconditioning techniques
            with dropping can be used. For the graphs related to some of the
            sparse matrices downloaded from the UF Sparse Matrix Collection, the
            experiments are compared to the known novel schemes. The results show
            that the provided method is more accurate. While it is slower than
            MC73 sequentially, it has good parallel efficiency compared with
            TraceMin.  In addition, it is insensitive to the selection of
            parameters, which is superior to the other two methods.
     ----------------------------------------------------------------------
                                                                     1,1
    Example:
    J=sorted(list(gb.Matrix.ssget(2760)), key=itemgetter(0))[0]
    J[1].print(level=3)
    pattern=J[1].pattern(gb.types.FP64)
    h=pattern+pattern.transpose()
    omega=lap.laplacian(pattern)
    #omega.print(level=3)
    hdip2=hdip_fiedler(omega)
    print(hdip2)

    The first line is simply the way to get a matrix from the SuiteSparse collection in python using the
    id of a matrix. It requires the import of "itemgetter" from operator, and the prior hand installation
    of ssgetpy.

    hdip2 from the example contains the vector as the first element in the tuple, the lambda value as the second element
    with the number of iterations as a tuple with a tuple in the third element.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
    """

    # Set u(0) = 1+sqrt(n),  u(1:n) = 1, and alpha = n+sqrt(n)
    n=L[0].nrows
    u=gb.Vector.dense(gb.types.FP64,n,fill=1)
    u[0]=1+sqrt(n)
    alpha=n+sqrt(n)

    #Set x(0) = 0 and x(1:n) = 1
    x=gb.Vector.dense(gb.types.FP64,n,fill=1)
    x[0]=0

    #set lnor equal to the inf-norm of L:
    lnor=L[1]

    #indiag = diagonal matrix with c = 1/L[0], as a preconditioner
    k_inner=0
    c=L[0].diag()
    indiag=(1/c.cast(gb.FP64)).diag()
    last_err=float('inf')
    #for i from 1 to kmax+1
    for i in range(1,kmax[0]+1):
        #compute beta = 2-norm of x  and set x = x/beta
 	x[0]=0
        #norm2, is a function to compute 2-norm of a vector imported from python file norm2, created by the same authors.
        beta = norm2(x)
        x=x/beta
        #using hmhx, a function to compute y = H*M*H*x
        y=hmhx(L[0],u,x,alpha)
        y[0]=0
        lamb= x.emult(y).reduce_float()
        print(lamb)
        #getting the inf norm for the vector normer, using norm(v,inf) = max(sum(abs(v))) vector v
        normer=y-lamb*x
        e=(normer.apply(gb.types.FP64.ABS)).max()
        # if e / inf norm of L < emax, break
        k_out=i
        #lnor is correct
        err=e/lnor
        print("i: "+str(i)+" err: "+str(err))
        if(err<emax or last_err<2*err):
            print("broken")
            break
        last_err=err
        #x(1:n) = L2 \ x(1:n)
        my=mypcg2(L[0],u,alpha,indiag,x,tol,kmax[1])
        x=my[0]
        kk=my[1]
        k_inner = k_inner+kk
        x[0]=0
    #compute beta = u'*x/alpha
    beta = (u.emult(x).reduce_float())/alpha
    #set x = x - beta*u
    x=x-beta*u
    #iters is used to return the number of inner and outer iterations
    iters=[k_out,k_inner]
    return [x,lamb,iters]
