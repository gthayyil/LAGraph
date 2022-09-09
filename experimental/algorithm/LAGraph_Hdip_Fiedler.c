#include "LG_internal.h"
#include <float.h>

//------------------------------------------------------------------------------
// LAGraph_Hdip_Fiedler
//------------------------------------------------------------------------------

//TODO:work on freeing the variables after implementation of hdip_fiedler
//-------------------------------------------------------------------------------------------------
int LAGraph_Happly //happly
(
    //outputs:
    GrB_Vector y, // y output of Householder reflection on x.
    //inputs:
    GrB_Vector u, // u, the vector used for application of householder
    GrB_Vector x, // x, the vector on which householder reflection is applied
    float alpha // the scalar alpha used for application of householder
)
{
    float reduced ;
    GRB_TRY (GrB_eWiseAdd (y,NULL, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32, u, x, NULL));
    GRB_TRY (GrB_reduce(&reduced, NULL,GrB_PLUS_FP32,y,NULL));
    GRB_TRY (GrB_apply(y,NULL,NULL,GrB_TIMES_FP32,-reduced/alpha,y,NULL));
    GRB_TRY (GrB_eWiseAdd(y,NULL,NULL,GrB_PLUS_FP32,GrB_PLUS_FP32,x,y,NULL));
    return (GrB_SUCCESS) ;
} 
//-------------------------------------------------------------------------------------------------
int LAGraph_hmhx
(
    //outputs:
    GrB_Vector z, //z output of hmhx
    //inputs:
    GrB_Matrix M, //Matrix used in hmhx
    GrB_Vector u, // Vector u used for happly
    GrB_Vector x, // Vector x used for happly 
    float alpha // the scalar alpha used for happly
) 
{
    GRB_TRY (LAGraph_Happly(z,u,x,alpha));
    GRB_TRY (GrB_mxv(z,NULL,NULL,NULL,M,z,NULL)); // matrix vector multiply with z and M
    GRB_TRY (LAGraph_Happly(z,u,z,alpha));// z  = happly with z,u and alpha
    GRB_TRY (GrB_Vector_setElement_FP32(z, 0, 0));
    return (GrB_SUCCESS) ;
}
//-------------------------------------------------------------------------------------------------
int LAGraph_norm2
(
    float *norm2,
    GrB_Vector v
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
    GRB_TRY (GrB_reduce (norm2, NULL, NULL, GrB_PLUS_FP32, t, NULL)) :
    &norm2 = sqrtf (&norm2) ;
    GrB_free (&t) ;
    return (GrB_SUCCESS) ;
}

//-------------------------------------------------------------------------------------------------
int LAGraph_Laplacian   // compute the Laplacian matrix of G->A
(
    // outputs:
    GrB_Matrix *Lap,    // Laplacian of G->A
    float *inform,      // infinity norm of Lap
    // inputs:
    LAGraph_Graph G,    // input matrix, symmetric
//  GrB_Type type,      // the type of Lap, typically GrB_FP32, ...
    char *msg
)
{
    GrB_Index ncol;
    GrB_Vector k = NULL ;
    GrB_Vector t = NULL ;
    GrB_Matrix *sparseM;
    GrB_Matrix *DMatrix;
    float *inform;
    
    // Assert G->A is symmetric
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->structure_is_symmetric == LAGraph_TRUE))
    {
    // Lap = (float) offdiag (G->A)
    GRB_TRY (GrB_Matrix_ncols(&ncol, G->A)) ;
    GRB_TRY (GrB_Matrix_new (Lap, GrB_FP32, ncol, ncol)) ;
    GRB_TRY (GrB_select (*Lap,NULL,NULL,GrB_OFFDIAG,G->A,0,NULL,NULL)) ;
    }
    else
    {
        // A is not known to be symmetric
        LG_ASSERT_MSG (false, -1005, "G->A must be symmetric") ;
    }

    // no self edges can be present
    LG_ASSERT_MSG (G->ndiag == 0, -1004, "G->ndiag must be zero") ;

    // t = row degree of Lap

    // t = Lap * x via the LAGraph_plus_one_fp32 semiring
    GRB_TRY (GrB_Vector_new (&t, GrB_FP32, ncol)) ;
    GRB_TRY(GrB_mxv(t,NULL,GrB_FP32,LAGraph_plus_one_fp32,*Lap,t,NULL));


    //creates a sparse Matrix with same dimensions as &Lap, and assigns -1 with &Lap as a Mask
    GRB_TRY (GrB_Matrix_new (sparseM, GrB_FP32, ncol, ncol)) ;
    //Python code has descriptor = S indicating structural mask
    GRB_TRY(GrB_assign(*sparseM,*Lap,NULL,-1,NULL,NULL,NULL,NULL, GrB_DESC_SC);

    
    //create a mask of 0s in vector t, and use that to replace the 0s with 1s. 
    GRB_TRY (GrB_Vector_new (&k, GrB_FP32, ncol)) ;
    GRB_TRY(GxB_select(k,NULL,NULL,GxB_EQ_ZERO,t,NULL,NULL);
    //Python code has descriptor = S indicating structural mask
    GRB_TRY(GrB_assign(t,k,NULL,1,NULL,NULL,GrB_DESC_SC);

    //inf norm calc using vector d and MAX_MONOID 
    GRB_TRY (GrB_reduce (inform, NULL, GrB_MAX_MONOID, t, NULL));
    *inform=*inform*2;
    
    //Using Matrix_diag to create a diagonal matrix from a vector    
    GRB_TRY (GrB_Matrix_new (DMatrix, GrB_FP32, ncol, ncol)) ;
    GRB_TRY (GrB_Matrix_diag(*DMatrix,t,NULL,NULL));
    
    //Calculating the Laplacian by adding the Dmatrix with SparseM.    
    GRB_TRY (GrB_eWiseAdd (*Lap, NULL, NULL, NULL, *DMatrix, *sparseM, NULL)) :
    //FREE k and sparseM

    return (GrB_SUCCESS);

int LAGraph_Hdip_Fiedler   // compute the Hdip_Fiedler
(
    // outputs:
    GrB_Vector iters, //This is a vector with number of inner and outer iterations listed in that order.
    float *lamb,    // Lambda of hdip_fiedler
    float *x,      // infinity norm of Lap
    // inputs:
    LAGraph_Matrix *L,    // input matrix, symmetric, result from Laplacian
    float *InfNorm,
    LAGraph_Vector kmax,
    float *emax,
    float *tol,
//  GrB_Type type,      // the type of Lap, typically GrB_FP32, ...
    char *msg
)
{   
    //TODO: Setup default values for emax, tol, and kmax
    //TODO: Setup freeing memory.

    GrB_Index n;
    GrB_Vector u = NULL ;
    GrB_Vector x = NULL ;
    GrB_Vector y = NULL ;
    GrB_Vector lambhelper = NULL;
    GrB_Matrix *c; // used to store diagonal of L
    float *alpha;
    float *k_inner;
    float *k_outer;
    float *last_err;
    float *beta;
    float *e;
    int i; // This is the integer used in for loop
    int kmaxZero; // kmax[0]

    //Set u(0) = 1+sqrt(n). u(1:n) = 1 and alpha = n+sqrt(n)
    GRB_TRY (GrB_Matrix_nrows(&n, *L)) ;
    GRB_TRY (GrB_Vector_new (&u, GrB_FP32, n)) ;
    GRB_TRY (GrB_apply (u, NULL, NULL, GxB_ONE_FP32, NULL, NULL)) ;
    GRB_TRY (GrB_Vector_setElement_FP32(u, 1+sqrtf(n), 0)) ;
    *alpha = n+sqrtf(n);

    //Set x(0) = 0 and x(1:n) = 1
    GRB_TRY (GrB_Vector_new (&x, GrB_FP32, n)) ;
    GRB_TRY (GrB_apply (x, NULL, NULL, GxB_ONE_FP32, NULL, NULL)) ;
    GRB_TRY (GrB_Vector_setElement_FP32(x, 0, 0)) ;

    //indiag = diagonal matrix with c = 1/L[0], as a preconditioner
    
    GRB_TRY (GrB_Matrix_new(c, GrB_FP32, n, n)) ;
    GRB_TRY (GrB_select (*c,NULL,NULL,GrB_DIAG,NULL,NULL,NULL,NULL)) ;
    GRB_TRY (GrB_apply (*c, NULL, NULL, GrB_DIV_FP32, *c,1, NULL)) ;
    *last_err = FLT_MAX;    

    //for i from 1 to kmax+1
    GRB_TRY(GrB_Vector_extractElement(&kmaxZero,kmax,0));
    for (i=1;i<=kmaxZero;k++)
    {
        //compute beta = 2-norm of x and set x = x/beta
        GRB_TRY (GrB_Vector_setElement_FP32(x, 0, 0)) ;
        GRB_TRY (LAGraph_norm2(beta,x));// z  = happly with z,u and alpha
        GRB_TRY (GrB_apply (x, NULL, NULL, GrB_RDIV_FP32, *beta,x, NULL)) ;
                
        //Set y = hmhx with m being L
        GRB_TRY (GrB_Vector_new (&y, GrB_FP32, n));
        GRB_TRY (LAGraph_hmhx(y,*L,u,x,alpha)); 
        GRB_TRY (GrB_Vector_setElement_FP32(y, 0, 0));
        //lamb = x.emult(y).reduce_float(mon=gb.types.FP32.PLUS_MONOID)
        GRB_TRY (GrB_Vector_new (&lambhelper, GrB_FP32, n));
        GRB_TRY (GrB_eWiseMult (lambhelper, NULL, NULL, GrB_TIMES_FP32, x, y, NULL)) :
        GRB_TRY (GrB_reduce(lamb,NULL,GrB_PLUS_FP32,lambhelper,NULL));
        GRB_TRY (GrB_Vector_clear(lambhelper));

        //getting the inf norm for the vector normer using norm(v,inf) = max(sum(abs(v))) vector v
        GRB_TRY (GrB_apply (lambhelper, NULL, NULL, GrB_TIMES_FP32,-lamb,x, NULL)) ;
        GRB_TRY (GrB_eWiseAdd (lambhelper, NULL, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32,y,lambhelper,NULL);
        //getting abs(lambhelper)
        GRB_TRY (GrB_apply (lambhelper, NULL, NULL, GrB_ABS_FP32,lambhelper, NULL)) ;
        GRB_TRY (GrB_apply (lambhelper, NULL, NULL, GrB_MAX_FP32,lambhelper, NULL)) ;
        GRB_TRY(GrB_Vector_extractElement(e,lambhelper,0));
        
        //if e/inf norm of L<emax, break
        *k_outer = i;
        *e = *e / *InfNorm
        
        if(*e<emax && *last_err<2*(*e))
        {
            break;
        }
        *lasterr=*e;
        //TODO: IMPLEMENT MYPCG2
    }
    
    //lamb = u.emult(x).reduce_float(mon=gb.types.FP32.PLUS_MONOID)
    GRB_TRY (GrB_Vector_clear(lambhelper));
    GRB_TRY (GrB_eWiseMult (lambhelper, NULL, NULL, GrB_TIMES_FP32, u, x, NULL)) :
    GRB_TRY (GrB_reduce(beta,NULL,GrB_PLUS_FP32,lambhelper,NULL));
    *beta = (*beta)/(*alpha);
       
       
    GRB_TRY (GrB_Vector_clear(lambhelper));
    GRB_TRY (GrB_apply (lambhelper, NULL, NULL, GrB_TIMES_FP32,-beta,u, NULL)) ;
    GRB_TRY (GrB_eWiseAdd (x, NULL, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32,x,lambhelper,NULL);

    //TODO: Do vectors start at 0 or 1?
    GRB_TRY (GrB_Vector_new (&iters, GrB_FP32, 2));
    GRB_TRY (GrB_Vector_setElement(iters,k_inner,0));
    GRB_TRY (GrB_Vector_setElement(iters,k_outer,1));

    return (GrB_SUCCESS);
}



