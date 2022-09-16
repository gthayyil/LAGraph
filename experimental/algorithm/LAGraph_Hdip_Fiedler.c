#include "LG_internal.h"
#include <float.h>

//------------------------------------------------------------------------------
// LAGraph_Hdip_Fiedler
//------------------------------------------------------------------------------

//TODO:work on freeing the variables after implementation of hdip_fiedler
//-------------------------------------------------------------------------------------------------
/*
The happly function, applies a Householder Reflection
    y= happly(u,x,alpha)
    calculates y = H*x where H = (I - u*u'/alpha)

    H*x =(I - u*u'/alpha)*x=(x-u*(u'u*x)/alpha)

    y,u, and x are vectors of size n and alpha is a scalar.
    The inputs are all modified so y must be a different vector.

    Example:
    I=[4.0,5.0,6.0,7.0]
    v=gb.Vector.from_list(I)
    J=[1.0,2.0,3.0,4.0]
    x=gb.Vector.from_list(J)
    A = 2.0
    k=happly(v,x,A)
    print(k)

    This example, simply creates 2 vectors and uses a literal numerical value as alpha to calculate happly.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm

*/
int LAGraph_Happly //happly Checked for pointer issues
(
    //outputs:
    GrB_Vector y, // y output of Householder reflection on x.
    //inputs:
    GrB_Vector u, // u, the vector used for application of householder
    GrB_Vector x, // x, the vector on which householder reflection is applied
    float alpha, // the scalar alpha used for application of householder
    //error msg
    char *msg
)
{
    float reduced;
    GRB_TRY (GrB_eWiseAdd (y,NULL, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32, u, x, NULL));
    GRB_TRY (GrB_reduce(reduced, NULL,GrB_PLUS_FP32,y,NULL));
    GRB_TRY (GrB_apply(y,NULL,NULL,GrB_TIMES_FP32,-reduced/alpha,y,NULL));
    GRB_TRY (GrB_eWiseAdd(y,NULL,NULL,GrB_PLUS_FP32,GrB_PLUS_FP32,x,y,NULL));
    
    GrB_free(&reduced);
    return (GrB_SUCCESS) ;
} 
//-------------------------------------------------------------------------------------------------
/*
hmhx = a function used to compute y = H*M*H*x = (M-u*x'-x*u)*x where x and y have the same size n

    Example:
    I=[4.0,5.0,6.0,7.0]
    v=gb.Vector.from_list(I)
    x=gb.Vector.from_list(I)
    A = 2.0
    I=[3,3,3]
    J=[1,2,3]
    K=[4,5,6]
    M=gb.Matrix.from_lists(I,J,K)
    M=M+M.transpose()
    mdiag=M.diag()
    indiag=1/mdiag.cast(gb.FP64)
    print(hmhx(indiag,v,x,A))

    This example starts off by creating 2 vectors and a matrix.
    Then makes sure the matrix is symmetric, and  gets the inverse diagonal of the matrix.
    Lastly, the example prints, hmhx called using the modified matrix, and vectors.


    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm

*/
int LAGraph_hmhx //hmhx checked for pointer issues
(
    //outputs:
    GrB_Vector z, //z output of hmhx
    //inputs:
    GrB_Matrix M, //Matrix used in hmhx
    GrB_Vector u, // Vector u used for happly
    GrB_Vector x, // Vector x used for happly 
    float alpha, // the scalar alpha used for happly
    char *msg
) 
{
    LG_TRY (LAGraph_Happly(z,u,x,alpha,msg));
    GRB_TRY (GrB_mxv(z,NULL,NULL,NULL,M,z,NULL)); // matrix vector multiply with z and M
    LG_TRY (LAGraph_Happly(z,u,z,alpha,msg));// z  = happly with z,u and alpha
    GRB_TRY (GrB_Vector_setElement_FP32(z, 0, 0));
    
    return (GrB_SUCCESS) ;
}
//-------------------------------------------------------------------------------------------------
/*
   norm2: This function aims to get the 2 norm of a vetor, and requires an import from python's math module.

    Example:
    I=[4.0,5.0,6.0,7.0]
    v=gb.Vector.from_list(I)
    k=norm2(v)
    print(k)

    Creates a vector, and passes it to the function and gets the norm2 of the function as a result.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
*/
int LAGraph_norm2 //norm2 checked for pointer mistakes
(
    //outputs:
    float norm2,
    //inputs:
    GrB_Vector v,
    //error msg
    char *msg
    
)
{
    GrB_Vector t = NULL ;
    GrB_Index len ;

    GRB_TRY (GrB_Vector_size (&len, v)) ;
    GRB_TRY (GrB_Vector_new (&t, GrB_FP32, len)) ;
#if LG_SUITESPARSE
    // t = v.^2
    GRB_TRY (GrB_apply (t, NULL, NULL, GxB_POW_FP32, v, (float) 2, NULL)) ;
#else
    // t = t.*t
    GRB_TRY (GrB_eWiseMult (t, NULL, NULL, GrB_TIMES_FP32, t, t, NULL)) :
#endif
    GRB_TRY (GrB_reduce (&norm2, NULL, NULL, GrB_PLUS_FP32, t, NULL)) :
    norm2 = sqrtf (norm2) ;

    GrB_free (&t);
    GrB_free (&len);
    return (GrB_SUCCESS) ;
}

//-------------------------------------------------------------------------------------------------
/*
  Laplacian: This function computes the laplacian of a Matrix (This matrix has to be symmetric (a normal matrix, when added with its transpose becomes symmetric)), and returns a tuple with
    both the laplacian of the matrix, and the infinity norm of the matrix.

    Example:
    J=sorted(list(gb.Matrix.ssget(2760)), key=itemgetter(0))[0]
    h=J[1]+J[1].transpose()
    L=laplacian(h)
    print(L)

    The first line is simply the way to get a matrix from the SuiteSparse collection in python using the
    id of a matrix. It requires the import of "itemgetter" from operator, and the prior hand installation
    of ssgetpy.

    The second line simply makes sure that the matrix is symmetric by adding it with its transpose.

    The third line calls this laplacian function on the matrix and returns a tuple with the laplacian and the inf norm.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
*/
int LAGraph_Laplacian   // compute the Laplacian matrix  
(
    //Question/TODO:Should i change matricies from pointers to normal variable names?
    // outputs:
    GrB_Matrix *Lap,    // the output Laplacian matrix 
    float inform,      // infinity norm of Lap
    // inputs:
    GrB_Matrix *G,    // input matrix, symmetric
//  GrB_Type type,      // the type of Lap, typically GrB_FP32, ...
    char *msg
)
{
    GrB_Index ncol;
    GrB_Vector k = NULL;
    GrB_Vector t = NULL;
    GrB_Matrix *sparseM;
    GrB_Matrix *DMatrix;
    
    // Assert G->A is symmetric
    // Lap = (float) offdiag (G->A)
    GRB_TRY (GrB_Matrix_ncols(&ncol, *G));
    GRB_TRY (GrB_Matrix_new (Lap, GrB_FP32, ncol, ncol));
    GRB_TRY (GrB_select (*Lap,NULL,NULL,GrB_OFFDIAG,G->A,0,NULL,NULL));

    // TODO: ASSERT NO SELF EDGES

    // t = row degree of Lap

    // t = Lap * x via the LAGraph_plus_one_fp32 semiring
    GRB_TRY (GrB_Vector_new (&t, GrB_FP32, ncol));
    GRB_TRY(GrB_mxv(t,NULL,GrB_FP32,LAGraph_plus_one_fp32,*Lap,t,NULL));


    //creates a sparse Matrix with same dimensions as &Lap, and assigns -1 with &Lap as a Mask
    GRB_TRY (GrB_Matrix_new (sparseM, GrB_FP32, ncol, ncol));
    //Python code has descriptor = S indicating structural mask
    GRB_TRY(GrB_assign(*sparseM,*Lap,NULL,-1,NULL,NULL,NULL,NULL, GrB_DESC_SC);

    
    //create a mask of 0s in vector t, and use that to replace the 0s with 1s. 
    GRB_TRY (GrB_Vector_new (&k, GrB_FP32, ncol));
    GRB_TRY(GxB_select(k,NULL,NULL,GxB_EQ_ZERO,t,NULL,NULL);
    //Python code has descriptor = S indicating structural mask
    GRB_TRY(GrB_assign(t,k,NULL,1,NULL,NULL,GrB_DESC_SC);

    //inf norm calc using vector d and MAX_MONOID 
    GRB_TRY (GrB_reduce (*inform, NULL, GrB_MAX_MONOID, t, NULL));
    inform=inform*2;
    
    //Using Matrix_diag to create a diagonal matrix from a vector    
    GRB_TRY (GrB_Matrix_new (DMatrix, GrB_FP32, ncol, ncol));
    GRB_TRY (GrB_Matrix_diag(*DMatrix,t,NULL,NULL));
    
    //Calculating the Laplacian by adding the Dmatrix with SparseM.    
    GRB_TRY (GrB_eWiseAdd (*Lap, NULL, NULL, NULL, *DMatrix, *sparseM, NULL));
    
    //FREE everything
    GrB_free (&ncol);
    GrB_free (&k);
    GrB_free (&t);
    GrB_free (sparseM);
    GrB_free (DMatrix);
    
    return (GrB_SUCCESS);

//-------------------------------------------------------------------------------------------------
/*
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
*/
int LAGraph_mypcg2
(
    //outputs
    GrB_Vector steper,
    float k
    // inputs:
    GrB_Matrix *L,    // input matrix, symmetric, result from Laplacian
    GrB_Vector u, //vector u will be passed into another function to create Householder reflection
    float malpha,  //This float 
    GrB_Matrix *invdiag,
    GrB_Vector b,
    float tol,
    float maxit, 
    //error msging
    char *msg
)
{
  
    GrB_Vector r, // This vector will be a copy of vector b to make sure vector b remains unchanged.
    GrB_Vector z, //used to apply preconditioner 
    GrB_Vector rho_helper, // used to help calculate rho
    GrB_Vector p, //search direction
    GrB_Vector q, //used for hmhx after finding next step in direction p
    GrB_Vector gamma_helper,
    float n;
    float bsize;
    float rho;
    float rho_prior;
    float gamma;
    float stepsize;
    float rnorm;

    //Set n to be number of rows in Laplacian matrix
    GRB_TRY (GrB_Matrix_nrows(&n, *L));
    
    //Set bsize to number of entries in vector b    
    GRB_TRY (GrB_Vector_size(&bsize, b));
    
    //b[0] = 0 is required for the input
    GRB_TRY (GrB_Vector_setElement_FP32(b, 0, 0));
    
    //Set r to be equal to vector b
    GRB_TRY (GrB_Vector_new(&r, GrB_FP32, bsize));
    GRB_TRY (GrB_Vector_dup(&r, b));

    //Set steper to size n, and filled with 0s.
    GRB_TRY (GrB_Vector_new (&steper, GrB_FP32, n));
    GRB_TRY (GrB_apply (steper, NULL, NULL, GxB_FIRST_FP32, 0, NULL,NULL));
    
    //Initial rho
    rho=1;

    //-------------''' Definition for helper vectors '''---------------------------
    //define rho_helper   
    GRB_TRY (GrB_Vector_new (&rho_helper, GrB_FP32, bsize));
 
    //define the vector used to apply preconditioner
    GRB_TRY (GrB_Vector_new(&z, GrB_FP32, bsize));
    
    //define new search direction vector  
    GRB_TRY (GrB_Vector_new (&p, GrB_FP32, bsize));
    
    //define vector q used to do hmhx after next step in direction p
    GRB_TRY (GrB_Vector_new (&q, GrB_FP32, bsize));

    //define vector gamma_helper
    GRB_TRY (GrB_Vector_new (&gamma_helper, GrB_FP32, bsize));

    //------------------------------------------------------------------------------
    
    for (k=1;k<=maxit;k++)
    {
        //Apply the preconditioner, using hmhx
        LG_TRY (LAGraph_hmhx(z,*invdiag,u,r,malpha,msg)); 
        GRB_TRY (GrB_Vector_setElement_FP32(z, 0, 0));
     
        //save the prior rho
        rho_prior=rho;
        
        //rho = r.emult(z).reduce_float(mon=gb.types.FP32.PLUS_MONOID) compute new rho
        GRB_TRY (GrB_eWiseMult (rho_helper, NULL, NULL, GrB_TIMES_FP32, r, z, NULL));
        GRB_TRY (GrB_reduce(rho,NULL,GrB_PLUS_FP32,rho_helper,NULL));
        GRB_TRY (GrB_Vector_clear(rho_helper));

        if(k==1){
            //first step is the direction p=z
            GRB_TRY (GrB_Vector_dup(&p, z));
        }else{
            //subsequent step in direction p
            beta = rho/rho_prior;
            //p=beta*p
            GRB_TRY (GrB_apply (p, NULL, NULL, GrB_TIMES_FP32, p, beta, NULL)) ;
            //p=p+z
            GRB_TRY (GrB_apply (p, NULL, NULL, GrB_PLUS_FP32, p, z, NULL)) ;
        }
        //hmhx is used on q
        LG_TRY (LAGraph_hmhx(q,*L,u,p,malpha,msg)); 
        GRB_TRY (GrB_Vector_setElement_FP32(q, 0, 0));

        //gamma=p.emult(q).reduce_float(mon=gb.types.FP32.PLUS_MONOID)
        GRB_TRY (GrB_eWiseMult (gamma_helper, NULL, NULL, GrB_TIMES_FP32, p, q, NULL));
        GRB_TRY (GrB_reduce(gamma,NULL,GrB_PLUS_FP32,gamma_helper,NULL));
        GRB_TRY (GrB_Vector_clear(gamma_helper));
        
        //stepsize to take
        stepsize = rho/gamma;
        
        //takes a step towards the solution , steper = steper + stepsize*p
        GRB_TRY (GrB_apply (p, NULL, NULL, GrB_TIMES_FP32, p, stepsize, NULL));
        GRB_TRY (GrB_apply (steper, NULL, NULL, GrB_PLUS_FP32, steper, p, NULL));

        //r=r-stepsize*q 
        GRB_TRY (GrB_apply (r, NULL, NULL, GrB_MINUS_FP32, r, q, NULL));
        GRB_TRY (GrB_Vector_setElement_FP32(steper, 0, 0));
        GRB_TRY (GrB_Vector_setElement_FP32(r, 0, 0));

        
        LG_TRY (LAGraph_norm2(rnorm,r,msg));// z  = happly with z,u and alpha
        if(rnorm < tol){
            break;
        }   
    }

    return (GrB_SUCCESS);
}


//-------------------------------------------------------------------------------------------------
/*
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

    Example:
    from operator import itemgetter
    J=sorted(list(gb.Matrix.ssget(2760)), key=itemgetter(0))[0]
    J[1].print(level=3)
    pattern=J[1].pattern(gb.types.FP32)
    h=pattern+pattern.transpose()
    omega=lap.laplacian(pattern)
    #omega.print(level=3)
    hdip2=hdip_fiedler(omega[0],omega[1])
    print(hdip2)

    The first line is simply the way to get a matrix from the SuiteSparse collection in python using the
    id of a matrix. It requires the import of "itemgetter" from operator, and the prior hand installation
    
    hdip2 from the example contains the vector as the first element in the tuple, the lambda value as the second element
    with the number of iterations as a tuple with a tuple in the third element.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
*/
int LAGraph_Hdip_Fiedler   // compute the Hdip_Fiedler
(
    // outputs:
    GrB_Vector iters, //This is a vector with number of inner and outer iterations listed in that order.
    float lamb,    // Lambda of hdip_fiedler
    GrB_Vector x; // the hdip fielder result vector
    // inputs:
    GrB_Matrix *L,    // input matrix, symmetric, result from Laplacian
    float InfNorm,
    GrB_Vector kmax,
    float emax,
    float tol,
//  GrB_Type type,      // the type of Lap, typically GrB_FP32, ...
    char *msg
)
{   
    //TODO: Setup default values for emax, tol, and kmax
    //TODO: Setup freeing memory.

    GrB_Index n;
    GrB_Vector u = NULL;
    GrB_Vector y = NULL;
    GrB_Vector lambhelper = NULL;
    GrB_Matrix *c; // used to store diagonal of L
    float alpha;
    float k_inner;
    float k_outer;
    float last_err;
    float beta;
    float e;
    float kk; //used to hold output from mypcg2
    int i; // This is the integer used in for loop
    int kmaxZero; // kmax[0]
    int kmaxOne; // kmax[1]
    
    //Set u(0) = 1+sqrt(n). u(1:n) = 1 and alpha = n+sqrt(n)
    GRB_TRY (GrB_Matrix_nrows(&n, *L));
    GRB_TRY (GrB_Vector_new (&u, GrB_FP32, n));
    GRB_TRY (GrB_apply (u, NULL, NULL, GxB_ONE_FP32, u, NULL));
    GRB_TRY (GrB_Vector_setElement_FP32(u, 1+sqrtf(n), 0));
    alpha = n+sqrtf(n);

    //Set x(0) = 0 and x(1:n) = 1
    GRB_TRY (GrB_Vector_new (&x, GrB_FP32, n));
    GRB_TRY (GrB_apply (x, NULL, NULL, GxB_ONE_FP32, x, NULL));
    GRB_TRY (GrB_Vector_setElement_FP32(x, 0, 0));

    //indiag = diagonal matrix with indiag = 1/L[0], as a preconditioner
    GRB_TRY (GrB_Matrix_new(indiag, GrB_FP32, n, n));
    GRB_TRY (GrB_select (*indiag,NULL,NULL,GrB_DIAG,NULL,NULL,NULL,NULL));
    GRB_TRY (GrB_apply (*indiag, NULL, NULL, GrB_RDIV_FP32, *indiag,1, NULL));
    last_err = FLT_MAX;    

    //for i from 1 to kmax[0]+1
    GRB_TRY(GrB_Vector_extractElement(&kmaxZero,kmax,0));
    //setting up kmax[1]
    GRB_TRY(GrB_Vector_extractElement(&kmaxOne,kmax,1));
    for (i=1;i<=kmaxZero;k++)
    {
        //compute beta = 2-norm of x and set x = x/beta
        GRB_TRY (GrB_Vector_setElement_FP32(x, 0, 0));
        LG_TRY (LAGraph_norm2(beta,x,msg));// z  = happly with z,u and alpha
        GRB_TRY (GrB_apply (x, NULL, NULL, GrB_RDIV_FP32, beta,x, NULL));
                
        //Set y = hmhx with m being L
        GRB_TRY (GrB_Vector_new (&y, GrB_FP32, n));
        LG_TRY (LAGraph_hmhx(y,*L,u,x,alpha,msg)); 
        GRB_TRY (GrB_Vector_setElement_FP32(y, 0, 0));
        //lamb = x.emult(y).reduce_float(mon=gb.types.FP32.PLUS_MONOID)
        GRB_TRY (GrB_Vector_new (&lambhelper, GrB_FP32, n));
        GRB_TRY (GrB_eWiseMult (lambhelper, NULL, NULL, GrB_TIMES_FP32, x, y, NULL));
        GRB_TRY (GrB_reduce(lamb,NULL,GrB_PLUS_FP32,lambhelper,NULL));
        GRB_TRY (GrB_Vector_clear(lambhelper));

        //getting the inf norm for the vector normer using norm(v,inf) = max(sum(abs(v))) vector v
        GRB_TRY (GrB_apply (lambhelper, NULL, NULL, GrB_TIMES_FP32,-lamb,x, NULL)) ;
        GRB_TRY (GrB_eWiseAdd (lambhelper, NULL, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32,y,lambhelper,NULL);
        //getting abs(lambhelper)
        GRB_TRY (GrB_apply (lambhelper, NULL, NULL, GrB_ABS_FP32,lambhelper, NULL));
        GRB_TRY (GrB_apply (lambhelper, NULL, NULL, GrB_MAX_FP32,lambhelper, NULL));
        GRB_TRY(GrB_Vector_extractElement(e,lambhelper,0));
        
        //if e/inf norm of L<emax, break
        k_outer = i;
        e = e/InfNorm;
        
        if(e<emax && last_err<2*e)
        {
            break;
        }
        lasterr=e;
        //x=mypcg2(L,u,alpha,indiag,x,tol,kmax[1])
        LG_TRY (LAGRAPH_mypcg2(x,kk,*L,alpha,indiag,x,tol,kmaxOne,msg));
        k_inner=k_inner+kk
        GRB_TRY (GrB_Vector_setElement_FP32(x, 0, 0));
 
    }
    
    //beta = u.emult(x).reduce_float(mon=gb.types.FP32.PLUS_MONOID)
    GRB_TRY (GrB_Vector_clear(lambhelper));
    GRB_TRY (GrB_eWiseMult (lambhelper, NULL, NULL, GrB_TIMES_FP32, u, x, NULL));
    GRB_TRY (GrB_reduce(beta,NULL,GrB_PLUS_FP32,lambhelper,NULL));
    beta = beta/alpha;
       
   
    GRB_TRY (GrB_Vector_clear(lambhelper));
    GRB_TRY (GrB_apply (lambhelper, NULL, NULL, GrB_TIMES_FP32,-beta,u, NULL));
    GRB_TRY (GrB_eWiseAdd (x, NULL, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32,x,lambhelper,NULL);

    //vectors start at 0
    GRB_TRY (GrB_Vector_new (&iters, GrB_FP32, 2));
    GRB_TRY (GrB_Vector_setElement(iters,k_inner,0));
    GRB_TRY (GrB_Vector_setElement(iters,k_outer,1));

    return (GrB_SUCCESS);
}



