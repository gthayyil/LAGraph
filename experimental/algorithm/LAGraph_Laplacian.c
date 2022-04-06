
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

#define LAGraph_FREE_WORK           \
{                                   \
    GrB_free (&T) ;                 \
    GrB_free (&u) ;                 \
    GrB_free (&w) ;                 \
    GrB_free (&y) ;                 \
    GrB_free (&L) ;                 \
    GrB_free (&thunk) ;             \
}


    LAGraph_FREE_WORK ;             \
    GrB_free (centrality) ;         \
}

#include "LG_internal.h"

//------------------------------------------------------------------------------
// LAGraph_VertexCentrality_Triangle: vertex triangle-centrality
//------------------------------------------------------------------------------

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
int LAGraph_norm2
(
    float *norm2,
    GrB_Vector v
)
{
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
    *norm2 = sqrtf (*norm2) ;
    GrB_free (&t) ;
    return (GrB_SUCCESS) ;
}

int LAGraph_Laplacian       // vertex triangle-centrality
(
    // outputs:
    GrB_Matrix *Lap,     // centrality(i): triangle centrality of i
    float *inform ,       // # of triangles in the graph
    // inputs:
    LAGraph_Graph G,            // input matrix, symmetric
    char *msg
)
{

    GRB_TRY (GrB_Matrix_ncols(&ncol, Lap)) ;
    GRB_TRY (GrB_Vector_new (&t, GrB_FP32, ncol)) ;
    GRB_TRY(GxB_select(Lap,NULL,NULL,GxB_OFFDIAG,NULL,NULL,NULL))
    GRB_TRY(GrB_mxv(t,NULL,NULL,NULL,))
    return (GrB_SUCCESS);
}
