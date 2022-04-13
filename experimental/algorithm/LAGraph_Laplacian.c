
//------------------------------------------------------------------------------
// LAGraph_VertexCentrality_triangle: vertex triangle-centrality
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University, ...

//------------------------------------------------------------------------------

// LAGraph_VertexCentrality_Triangle: computes the TriangleCentrality of
// an undirected graph.  No self edges are allowed on the input graph.
// Methods 2 and 3 can tolerate any edge weights (they are ignored; only the
// pattern of G->A is used).  Methods 1 and 1.5 require unit edge weights
// (this could be modified); results are undefined if this condition doesn't
// hold.

// P. Burkhardt, "Triangle centrality," https://arxiv.org/pdf/2105.00110.pdf,
// April 2021.

// Methods 2 and 3 require SuiteSparse:GraphBLAS.  Method 3 is by far the
// fastest.

// TC0: in python (called TC1 in the first draft of the paper)
//
// def triangle_centrality1(A):
//          T = A.mxm(A, mask=A)
//          y = T.reduce_vector()
//          k = y.reduce_float()
//          return(1/k)*(3*(A @ y) - 2*(T @ y) + y)
//          note: T@y is wrong. should be plus_second semiring

//  def TC1(A):
//      # this was "Method 1.5" in a draft, note the T.one@y is now correct:
//      T = A.mxm(A, mask=A, desc=ST1)
//      y = T.reduce_vector()
//      k = y.reduce_float()
//      return (3 * (A @ y) - 2 * (T.one() @ y) + y) / k

//  def TC2(A):
//      # this was TC2 in the first submission
//      T = A.plus_pair(A, mask=A, desc=ST1)
//      y = Vector.dense(FP64, A.nrows)
//      T.reduce_vector(out=y, accum=FP64.plus)
//      k = y.reduce_float()
//      return (3 * A.plus_second(y) - 2 * T.plus_second(y) + y) / k

//  def TC3(A):
//      M = A.tril(-1)
//      T = A.plus_pair(A, mask=M, desc=ST1)
//      y = T.reduce() + T.reduce(desc=ST0)
//      k = y.reduce_float()
//      return (
//          3 * A.plus_second(y) -
//          (2 * (T.plus_second(y) + T.plus_second(y, desc=ST0))) + y
//      ) / k

//------------------------------------------------------------------------------

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&T) ;                 \
    GrB_free (&u) ;                 \
    GrB_free (&w) ;                 \
    GrB_free (&y) ;                 \
    GrB_free (&L) ;                 \
    GrB_free (&thunk) ;             \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (centrality) ;         \
}

#include "LG_internal.h"

//------------------------------------------------------------------------------
// LAGraph_Laplacian
//------------------------------------------------------------------------------

//TODO:work on freeing the variables after implementation of hdip_fiedler

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
    GrB_TRY (GrB_eWiseAdd (y, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32, u, x, NULL));
    GrB_TRY (GrB_reduce(&reduced, NULL,NULL,GrB_PLUS_FP32,y,NULL));
    GrB_TRY (GrB_apply(y,NULL,NULL,GrB_TIMES_FP32,-reduced/alpha,y,NULL));
    GrB_TRY (GrB_eWiseAdd(y,NULL,NULL,GrB_PLUS_FP32,GrB_PLUS_FP32,x,y,NULL));
    return 0;
} 

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
    GrB_TRY (LAGraph_Happly(z,u,x,alpha));
    GrB_TRY (GrB_mxv(z,NULL,NULL,NULL,M,z,NULL));
    GrB_TRY (LAGraph_Happly(z,u,z,alpha));
//    GrB_TRY (z[0]=0), use function from GrB
    GRB_TRY (GrB_Vector_setElement_FP32(z, 0, 0));
    return 0;
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL ;

int LAGraph_norm2
(
    float *norm2,
    GrB_Vector v,
    GrB_Vector t        // temporary vector, type GrB_FP32, same size as v
)
{
//  GrB_Vector t = NULL ;
//  GrB_Index len ;

//  GRB_TRY (GrB_Vector_size (&len, v)) ;
//  GRB_TRY (GrB_Vector_new (&t, GrB_FP32, len)) ;

#if LG_SUITESPARSE
    // t = v.^2
    GRB_TRY (GrB_apply (t, NULL, NULL, GxB_POW_FP32, v, (float) 2, NULL)) ;
#else
    // t = t.*t
    GRB_TRY (GrB_eWiseMult (t, NULL, NULL, GrB_TIMES_FP32, t, t, NULL)) :
#endif
    GRB_TRY (GrB_reduce (norm2, NULL, NULL, GrB_PLUS_FP32, t, NULL)) :
    *norm2 = sqrtf (*norm2) ;
//  GrB_free (&t) ;
    return (GrB_SUCCESS) ;
}

#undef  LG_FREE_WORK
#define LG_FREE_WORK                \
{                                   \
    GrB_free (&T) ;                 \
    GrB_free (&u) ;                 \
    GrB_free (&w) ;                 \
    GrB_free (&y) ;                 \
    GrB_free (&L) ;                 \
    GrB_free (&thunk) ;             \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (centrality) ;         \
}

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
    GrB_Matrix sparseM = NULL, DMatrix = NULL ;
    GrB_Vector t = NULL, k = NULL ;
    // TODO assert Lap and inform are not null
    (*Lap) = NULL ;

    // TODO: assert G->A is symmetric

    // Lap = (float) offdiag (G->A)
    GRB_TRY (GrB_Matrix_ncols(&ncol, G->A)) ;
    GRB_TRY (GrB_Matrix_new (Lap, GrB_FP32, ncol, ncol)) ;
    GRB_TRY (GrB_select (*Lap,NULL,NULL,GrB_OFFDIAG,G->A,0,NULL,NULL)) ;

    // t = row degree of Lap

    // t = Lap * x via the LAGraph_plus_one_fp32 semiring
    GRB_TRY (GrB_Vector_new (&t, GrB_FP32, ncol)) ;
    GRB_TRY(GrB_mxv(t,NULL,GrB_FP32,LAGraph_plus_one_fp32,*Lap,t,NULL));


    //creates a sparse Matrix with same dimensions as *Lap, and assigns -1 with *Lap as a Mask
    GRB_TRY (GrB_Matrix_new (&sparseM, GrB_FP32, ncol, ncol)) ;
    //Python code has descriptor = S, but im not sure its purpose
    GRB_TRY(GrB_assign(sparseM,*Lap,NULL,-1,NULL,NULL,NULL,NULL, GrB_DESC_S  );

    
    //create a mask of 0s in vector t, and use that to replace the 0s with 1s. 
    GRB_TRY (GrB_Vector_new (&k, GrB_FP32, ncol)) ;
    // TODO: use GrB_select instead
    GRB_TRY(GxB_select(k,NULL,NULL,GxB_EQ_ZERO,t,NULL,NULL);
    //Python code uses descriptor=gb.descriptor.S , unsure of its purpose
    GRB_TRY(GrB_assign(t,k,NULL,1,NULL,NULL, GrB_DESC_S);

    //inf norm calc using vector d and MAX_MONOID 
    GRB_TRY (GrB_reduce (inform, NULL, GrB_MAX_MONOID_FP32, t, NULL));
    *inform=(*inform)*2;
    
    //Using Matrix_diag to create a diagonal matrix from a vector    
// old:
//  GRB_TRY (GrB_Matrix_new (DMatrix, GrB_FP32, ncol, ncol)) ;
//  GRB_TRY (GrB_Matrix_diag(*DMatrix,t,NULL,NULL));
// new:
    GRB_TRY (GrB_Matrix_diag(DMatrix,t,NULL,NULL));
    
    //Calculating the Laplacian by adding the Dmatrix with SparseM.    
    GRB_TRY (GrB_eWiseAdd (*Lap, NULL, NULL, NULL, DMatrix, sparseM, NULL)) :
    LG_FREE_WORK ;
    return (GrB_SUCCESS);
}
