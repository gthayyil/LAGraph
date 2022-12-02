//------------------------------------------------------------------------------
// LAGraphX.h: include file for LAGraph experimental code
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

//------------------------------------------------------------------------------

#ifndef LAGRAPHX_H
#define LAGRAPHX_H

#include <GraphBLAS.h>
#include <LAGraph.h>

//==============================================================================
// Experimental methods: in experimental/algorithm and experimental/utility
//==============================================================================

// Do not rely on these in production.  These methods are still under
// development, and is intended only for illustration not benchmarking.  Do not
// use for benchmarking, without asking the authors.

//------------------------------------------------------------------------------
// LAGraph_Random_*: Random number generator
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC
int LAGraph_Random_Init
(
    char *msg
) ;
LAGRAPH_PUBLIC
int LAGraph_Random_Finalize
(
    char *msg
) ;

#if defined ( COVERAGE )
// for testing only
LAGRAPH_PUBLIC bool random_hack ;
#endif

LAGRAPH_PUBLIC
int LAGraph_Random_Seed     // construct a random seed vector
(
    // input/output
    GrB_Vector Seed,    // vector of random number seeds, normally GrB_UINT64
    // input
    uint64_t seed,      // scalar input seed
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_Random_Next     // advance to next random vector
(
    // input/output
    GrB_Vector Seed,
    char *msg
) ;

LAGRAPH_PUBLIC
GrB_Info LAGraph_Random_Matrix    // random matrix of any built-in type
(
    // output
    GrB_Matrix *A,      // A is constructed on output
    // input
    GrB_Type type,      // type of matrix to construct
    GrB_Index nrows,    // # of rows of A
    GrB_Index ncols,    // # of columns of A
    double density,     // density: build a sparse matrix with
                        // density*nrows*cols values if not INFINITY;
                        // build a dense matrix if INFINITY.
    uint64_t seed,      // random number seed
    char *msg
) ;

//****************************************************************************
// binary file I/O
//****************************************************************************

// The LAGraph *.lagraph file consists of an ASCII JSON header, followed by
// one or more serialized "blobs" created by GrB_Matrix_serialize (or
// GxB_Matrix_serialize if using SuiteSparse:GraphBLAS).  The file can only be
// read back into LAGraph when using the same GraphBLAS library used to create
// it.

// To create a binary file containing one or more GrB_Matrix objects, the user
// application must first open the file f, create the ascii JSON header with
// LAGraph_SWrite_Header*, and then write one or more binary serialized
// GrB_Matrix blobs from  using LAGraph_SWrite_Matrix.

// Example:

/*
    // serialize the matrices A (of type GrB_FP64) and B (of type GrB_BOOL)
    void *Ablob, *Bblob ;
    GrB_Index Ablob_size, Bblob_size ;
    GxB_Matrix_serialize (&Ablob, &Ablob_size, A, NULL) ;
    GxB_Matrix_serialize (&Bblob, &Bblob_size, B, NULL) ;

    // open the file and write the JSON header
    FILE *f = fopen ("mymatrices.lagraph", "w") ;
    LAGraph_SWrite_HeaderStart (f, "mystuff", msg) ;
    LAGraph_SWrite_HeaderItem (f, LAGraph_matrix_kind, "A", "double", 0,
        Ablob_size, msg) ;
    LAGraph_SWrite_HeaderItem (f, LAGraph_matrix_kind, "B", "bool", 0,
        Bblob_size, msg) ;
    LAGraph_SWrite_HeaderEnd (f, msg) ;

    // write the matrices in binary
    LAGraph_SWrite_Item (f, Ablob, Ablob_size, msg) ;
    LAGraph_SWrite_Item (f, Bblob, Bblob_size, msg) ;

    fclose (f) ;
*/

typedef enum
{
    LAGraph_unknown_kind = -1,  // unknown kind
    LAGraph_matrix_kind = 0,    // a serialized GrB_Matrix
    LAGraph_vector_kind = 1,    // a serialized GrB_Vector (SS:GrB only)
    LAGraph_text_kind = 2,      // text (char *), possibly compressed
}
LAGraph_Contents_kind ;

typedef struct
{
    // serialized matrix/vector, or pointer to text, and its size
    void *blob ;
    size_t blob_size ;

    // kind of item: matrix, vector, text, or unknown
    LAGraph_Contents_kind kind ;

    // if kind is text: compression used
    // -1: none, 0: default for library, 1000: LZ4, 200x: LZ4HC:x
    int compression ;

    // name of the object
    char name [LAGRAPH_MAX_NAME_LEN+4] ;

    // if kind is matrix or vector: type name
    char type_name [LAGRAPH_MAX_NAME_LEN+4] ;
}
LAGraph_Contents ;

LAGRAPH_PUBLIC
int LAGraph_SWrite_HeaderStart  // write the first part of the JSON header
(
    FILE *f,                    // file to write to
    const char *name,           // name of this collection of matrices
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_SWrite_HeaderItem   // write a single item to the JSON header
(
    // inputs:
    FILE *f,                    // file to write to
    LAGraph_Contents_kind kind, // matrix, vector, or text
    const char *name,           // name of the matrix/vector/text; matrices from
                                // sparse.tamu.edu use the form "Group/Name"
    const char *type,           // name of type of the matrix/vector
    int compression,            // text compression method
    GrB_Index blob_size,        // exact size of serialized blob for this item
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_SWrite_HeaderItem   // write a single item to the JSON header
(
    // inputs:
    FILE *f,                    // file to write to
    LAGraph_Contents_kind kind, // matrix, vector, or text
    const char *name,           // name of the matrix/vector/text; matrices from
                                // sparse.tamu.edu use the form "Group/Name"
    const char *type,           // name of type of the matrix/vector
    // todo: text not yet supported by LAGraph_SWrithe_HeaderItem
    int compression,            // text compression method
    GrB_Index blob_size,        // exact size of serialized blob for this item
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_SWrite_HeaderEnd    // write the end of the JSON header
(
    FILE *f,                    // file to write to
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_SWrite_Item  // write the serialized blob of a matrix/vector/text
(
    // input:
    FILE *f,                // file to write to
    const void *blob,       // serialized blob from G*B_Matrix_serialize
    GrB_Index blob_size,    // exact size of the serialized blob
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_SRead   // read a set of matrices from a *.lagraph file
(
    FILE *f,                        // file to read from
    // output
    char **collection,              // name of collection (allocated string)
    LAGraph_Contents **Contents,    // array contents of contents
    GrB_Index *ncontents,           // # of items in the Contents array
    char *msg
) ;

LAGRAPH_PUBLIC
void LAGraph_SFreeContents      // free the Contents returned by LAGraph_SRead
(
    // input/output
    LAGraph_Contents **Contents,    // array of size ncontents
    GrB_Index ncontents
) ;

LAGRAPH_PUBLIC
int LAGraph_SSaveSet            // save a set of matrices from a *.lagraph file
(
    // inputs:
    char *filename,             // name of file to write to
    GrB_Matrix *Set,            // array of GrB_Matrix of size nmatrices
    GrB_Index nmatrices,        // # of matrices to write to *.lagraph file
//  todo: handle vectors and text in LAGraph_SSaveSet
    char *collection,           // name of this collection of matrices
    char *msg
) ;

int LAGraph_SLoadSet            // load a set of matrices from a *.lagraph file
(
    // input:
    char *filename,             // name of file to read from
    // outputs:
    GrB_Matrix **Set_handle,        // array of GrB_Matrix of size nmatrices
    GrB_Index *nmatrices_handle,    // # of matrices loaded from *.lagraph file
//  todo: handle vectors and text in LAGraph_SLoadSet
//  GrB_Vector **Set_handle,        // array of GrB_Vector of size nvector
//  GrB_Index **nvectors_handle,    // # of vectors loaded from *.lagraph file
//  char **Text_handle,             // array of pointers to (char *) strings
//  GrB_Index **ntext_handle,       // # of texts loaded from *.lagraph file
    char **collection_handle,   // name of this collection of matrices
    char *msg
) ;

LAGRAPH_PUBLIC
void LAGraph_SFreeSet           // free a set of matrices
(
    // input/output
    GrB_Matrix **Set_handle,    // array of GrB_Matrix of size nmatrices
    GrB_Index nmatrices         // # of matrices in the set
) ;

//****************************************************************************
// Algorithms
//****************************************************************************

//****************************************************************************
/**
 * Given a symmetric graph A with no-self edges, compute all k-trusses of A.
 *
 * @param[out]  Cset    size n, output k-truss subgraphs.
 * @param[out]  kmax    smallest k where k-truss is empty
 * @param[out]  ntris   Array of size n (on input), ntris [k] is num triangles in k-truss
 * @param[out]  nedges  Array of size n (on input), nedges [k] is num edges in k-truss
 * @param[out]  nstepss Array of size n (on input), nstepss [k] is num steps for k-truss
 * @param[in]   G       input graph, A, not modified.  Must be undirected
 *                      or directed with symmetric structure, no self edges.
 * @param[in,out] msg   any error messages.
 *
 * @retval GrB_SUCCESS      if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER if kmax, ntris, nedges, nsteps is NULL
 */
LAGRAPH_PUBLIC
int LAGraph_AllKTruss   // compute all k-trusses of a graph
(
    // outputs
    GrB_Matrix *Cset,   // size n, output k-truss subgraphs
    int64_t *kmax,      // smallest k where k-truss is empty
    int64_t *ntris,     // size max(n,4), ntris [k] is #triangles in k-truss
    int64_t *nedges,    // size max(n,4), nedges [k] is #edges in k-truss
    int64_t *nstepss,   // size max(n,4), nstepss [k] is #steps for k-truss
    // input
    LAGraph_Graph G,    // input graph
    char *msg
) ;

//****************************************************************************
/**
 * Given an undirected graph G with no-self edges, LAGraph_KTruss finds the
 * k-truss subgraph of G.
 *
 * @param[out]  C       k-truss subgraph, of type GrB_UINT32
 * @param[in]   G       input graph, not modified
 * @param[in]   k       the truss to find
 * @param[in,out] msg   any error messages.
 *
 * @retval GrB_SUCCESS      if completed successfully (equal or not)
 * @retval GrB_NULL_POINTER if C or C_type is NULL
 * @return Any GraphBLAS errors that may have been encountered
 */
LAGRAPH_PUBLIC
int LAGraph_KTruss      // compute the k-truss of a graph
(
    // outputs:
    GrB_Matrix *C,      // output k-truss subgraph, C
    // inputs:
    LAGraph_Graph G,    // input graph
    uint32_t k,         // find the k-truss, where k >= 3
    char *msg
) ;

//****************************************************************************
// Connected components
//****************************************************************************

/**
 * Determine connected components in an undirected graph.
 *
 * @param[out] result    array of component identifiers for each vertex (allocated
 *                       by the algorithm, ownership returned to caller).
 * @param[in]  A         the graph (symmetric)
 * @param[in]  sanitize  If true, test to ensure A is symmetric
 * @param[in,out] msg    any error messages.
 *
 * @retval GrB_SUCCESS      if completed successfully
 * @retval GrB_NULL_POINTER if result is NULL
 */
LAGRAPH_PUBLIC
int LAGraph_cc_lacc (
    GrB_Vector *result,
    GrB_Matrix A,
    bool sanitize,
    char *msg
) ;

//****************************************************************************
// Bellman Ford variants
//****************************************************************************

/**
 * Bellman-Ford single source shortest paths, returning just the shortest path
 * lengths.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[in]   A            matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_basic
(
    GrB_Vector *pd_output,
    const GrB_Matrix A,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning just the shortest path
 * lengths.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[in]   A            matrix for the graph (optional-ish)
 * @param[in]   AT           transpose of A (optional-ish)
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output is NULL or both A and AT are NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_basic_pushpull
(
    GrB_Vector *pd_output,
    const GrB_Matrix A,
    const GrB_Matrix AT,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning just the shortest path
 * lengths.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[in]   AT           transposed adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output or AT is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_basic_mxv
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    const GrB_Matrix AT,        //transposed adjacency matrix for the graph
    const GrB_Index s           //given index of the source
) ;

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_full
(
    GrB_Vector *pd_output,
    GrB_Vector *ppi_output,
    GrB_Vector *ph_output,
    const GrB_Matrix A,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_full1
(
    GrB_Vector *pd_output,
    GrB_Vector *ppi_output,
    GrB_Vector *ph_output,
    const GrB_Matrix A,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_full1a
(
    GrB_Vector *pd_output,
    GrB_Vector *ppi_output,
    GrB_Vector *ph_output,
    const GrB_Matrix A,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or A is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_full2
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,      //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
) ;

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree.
 *
 * @param[out]  pd_output    the pointer to the vector of distance (created internally)
 * @param[out]  ppi_output   the pointer to the vector of parent (created internally)
 * @param[out]  ph_output    the pointer to the vector of hops (created internally)
 * @param[in]   AT           transpose of the adjacency matrix for the graph
 * @param[in]   s            index of the source
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd_output, ppi_output, ph_output, or AT is NULL
 * @retval GrB_INVALID_VALUE  if A is not square, s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_full_mxv
(
    GrB_Vector *pd_output,
    GrB_Vector *ppi_output,
    GrB_Vector *ph_output,
    const GrB_Matrix AT,
    const GrB_Index s
) ;

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree (integer weights).
 *
 * @param[out]  pd       pointer to distance vector d, d(k) = shortest distance
 *                       between s and k if k is reachable from s
 * @param[out]  ppi      pointer to parent index vector pi, pi(k) = parent of
 *                       node k in the shortest path tree
 * @param[in]   s        index of the source
 * @param[in]   n        number of nodes
 * @param[in]   nz       number of edges
 * @param[in]   I        row index vector (size n)
 * @param[in]   J        column index vector (size nz)
 * @param[in]   W        weight vector (size nz), W(i) = weight of edge (I(i),J(i))
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd, ppi, I, J, or W is NULL
 * @retval GrB_INVALID_VALUE  if s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails.
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_pure_c
(
    int32_t **pd,

    int64_t **ppi,

    const int64_t s,
    const int64_t n,
    const int64_t nz,
    const int64_t *I,
    const int64_t *J,
    const int32_t *W
) ;

/**
 * Bellman-Ford single source shortest paths, returning both the path lengths
 * and the shortest-path tree (double weights).
 *
 * @param[out]  pd       pointer to distance vector d, d(k) = shortest distance
 *                       between s and k if k is reachable from s
 * @param[out]  ppi      pointer to parent index vector pi, pi(k) = parent of
 *                       node k in the shortest path tree
 * @param[in]   s        index of the source
 * @param[in]   n        number of nodes
 * @param[in]   nz       number of edges
 * @param[in]   I        row index vector (size n)
 * @param[in]   J        column index vector (size nz)
 * @param[in]   W        weight vector (size nz), W(i) = weight of edge (I(i),J(i))
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If pd, ppi, I, J, or W is NULL
 * @retval GrB_INVALID_VALUE  if s is not a valid vertex index
 * @retval GrB_OUT_OF_MEMORY  if allocation fails.
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 *
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_BF_pure_c_double
(
    double **pd,

    int64_t **ppi,

    const int64_t s,
    const int64_t n,
    const int64_t nz,
    const int64_t *I,
    const int64_t *J,
    const double  *W
) ;

//****************************************************************************
/**
 * Community detection using label propagation algorithm
 *
 * @param[out]  CDLP_handle  community vector
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   symmetric    denote whether the matrix is symmetric
 * @param[in]   sanitize     if true, verify that A is binary
 * @param[in]   itermax      max number of iterations (0 computes nothing)
 * @param[out]  t            array of two doubles allocated by caller:
 *                           [0]=sanitize time, [1]=cdlp time in seconds
 * @param[in,out] msg        any error messages.
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NULL_POINTER   If t or CDLP_handle is NULL
 * @retval GrB_INVALID_OBJECT If A is not stored in CSR format
 * @retval GrB_OUT_OF_MEMORY  if allocation fails.
 * @retval GrB_NO_VALUE       if A has a negative weight cycle
 */
LAGRAPH_PUBLIC
int LAGraph_cdlp
(
    GrB_Vector *CDLP_handle,
    const GrB_Matrix A,
    bool symmetric,
    bool sanitize,
    int itermax,
    double *t,
    char *msg
) ;

//****************************************************************************
/**
 * Sparse deep neural network inference. Performs ReLU inference using input
 * feature vectors Y0.
 *
 * @param[out]  Yhandle      Y, created on output
 * @param[in]   W            W [0..nlayers-1], each nneurons-by-nneurons
 * @param[in]   Bias         Bias [0..nlayers-1], diagonal nneurons-by-nneurons
 * @param[in]   nlayers      number of layers
 * @param[in]   Y0           input features: nfeatures-by-nneurons
 *
 * @retval GrB_SUCCESS         if completed successfully
 * @retval GrB_NOT_IMPLEMENTED vanilla version has not been implemented yet
 * @retval GrB_NULL_POINTER    If Yhandle, W, Bias, or Y0 is NULL
 * @retval GrB_DOMAIN_MISMATCH if type of Y0 is not FP32 or FP64, or the types of
 *                             W or Bias arent the same as Y0
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_dnn
(
    // output
    GrB_Matrix *Yhandle,
    // input: not modified
    GrB_Matrix *W,
    GrB_Matrix *Bias,
    int nlayers,
    GrB_Matrix Y0
) ;

//****************************************************************************
/**
 * Compute all-pairs shortest paths using Floyd-Warshall method
 *
 * @param[in]   G       input graph, with edge weights
 * @param[out]  D       output graph, created on output
 * @param[out]  D_type  type of scalar stored in D (see source for explanation)
 *
 * @retval GrB_SUCCESS         if completed successfully
 * @retval GrB_NOT_IMPLEMENTED vanilla version has not been implemented yet
 * @retval GrB_NULL_POINTER    If D or D_type is NULL
 * @retval GrB_INVALID_VALUE   If G is not square
 */
LAGRAPH_PUBLIC
GrB_Info LAGraph_FW
(
    const GrB_Matrix G,
    GrB_Matrix *D,
    GrB_Type   *D_type
) ;

//****************************************************************************
/**
 * Compute the local clustering coefficient for all nodes in a graph.
 *
 * @param[out]  LCC_handle   output vector holding coefficients
 * @param[in]   A            adjacency matrix for the graph
 * @param[in]   symmetric    denote whether the matrix is symmetric
 * @param[in]   sanitize     if true, verify that A is binary
 * @param[out]  t            array of two doubles
 *                           [0]=sanitize time, [1]=lcc time in seconds
 * @param[in,out] msg        any error messages.
 *
 * @retval GrB_SUCCESS        if completed successfully
 * @retval GrB_NOT_IMPLEMENTED vanilla version has not been implemented yet
 * @retval GrB_NULL_POINTER   If LCC_handle or LCC_type is NULL
 * @retval GrB_INVALID_VALUE  If A is not stored in CSR format
 */
LAGRAPH_PUBLIC
int LAGraph_lcc            // compute lcc for all nodes in A
(
    GrB_Vector *LCC_handle,     // output vector
    const GrB_Matrix A,         // input matrix
    bool symmetric,             // if true, the matrix is symmetric
    bool sanitize,              // if true, ensure A is binary
    double t [2],               // t [0] = sanitize time, t [1] = lcc time,
                                // in seconds
    char *msg
) ;

//****************************************************************************

LAGRAPH_PUBLIC
int LAGraph_msf
(
    GrB_Matrix *result, // output: an unsymmetrical matrix, the spanning forest
    GrB_Matrix A,       // input matrix
    bool sanitize,      // if true, ensure A is symmetric
    char *msg
) ;

//****************************************************************************

LAGRAPH_PUBLIC
int LAGraph_scc (
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    char *msg
) ;

//****************************************************************************
LAGRAPH_PUBLIC
int LAGraph_VertexCentrality_Triangle       // vertex triangle-centrality
(
    // outputs:
    GrB_Vector *centrality,     // centrality(i): triangle centrality of i
    uint64_t *ntriangles,       // # of triangles in the graph
    // inputs:
    int method,                 // 0, 1, 2, or 3
    LAGraph_Graph G,            // input graph
    char *msg
) ;

//****************************************************************************
LAGRAPH_PUBLIC
int LAGraph_MaximalIndependentSet       // maximal independent set
(
    // outputs:
    GrB_Vector *mis,            // mis(i) = true if i is in the set
    // inputs:
    LAGraph_Graph G,            // input graph
    uint64_t seed,              // random number seed
    GrB_Vector ignore_node,     // if NULL, no nodes are ignored.  Otherwise
                                // ignore_node(i) = true if node i is to be
                                // ignored, and not treated as a candidate
                                // added to maximal independent set.
    char *msg
) ;

LAGRAPH_PUBLIC
int LG_CC_FastSV5           // SuiteSparse:GraphBLAS method, with GxB extensions
(
    // output
    GrB_Vector *component,  // output: array of component identifiers
    // inputs
    LAGraph_Graph G,        // input graph, modified then restored
    char *msg
) ;

//------------------------------------------------------------------------------
// kcore algorithms
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC
int LAGraph_KCore_All
(
    // outputs:
    GrB_Vector *decomp,     // kcore decomposition
    uint64_t *kmax,
    // inputs:
    LAGraph_Graph G,            // input graph
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_KCore
(
    // outputs:
    GrB_Vector *decomp,     // kcore decomposition
    // inputs:
    LAGraph_Graph G,        // input graph
    uint64_t k,             //k level to compare to
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_KCore_Decompose
(
    // outputs:
    GrB_Matrix *D,              // kcore decomposition
    // inputs:
    LAGraph_Graph G,            // input graph
    GrB_Vector decomp,         // input decomposition matrix
    uint64_t k,
    char *msg
) ;

//------------------------------------------------------------------------------
// counting graphlets
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC
int LAGraph_FastGraphletTransform
(
    // outputs:
    GrB_Matrix *F_net,  // 16-by-n matrix of graphlet counts
    // inputs:
    LAGraph_Graph G,
    bool compute_d_15,  // probably this makes most sense
    char *msg
) ;

LAGRAPH_PUBLIC
int LAGraph_SquareClustering
(
    // outputs:
    GrB_Vector *square_clustering,
    // inputs:
    LAGraph_Graph G,
    char *msg
) ;

//------------------------------------------------------------------------------
// HDIP_Fiedler
//------------------------------------------------------------------------------
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
LAGRAPH_PUBLIC
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
);
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
LAGRAPH_PUBLIC
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
);
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

LAGRAPH_PUBLIC
int LAGraph_norm2 //norm2 checked for pointer mistakes
(
    //outputs:
    float norm2,
    //inputs:
    GrB_Vector v,
    //error msg
    char *msg

);

//-------------------------------------------------------------------------------------------------
/*
  Laplacian: This function computes the laplacian of a Matrix (This matrix has to be symmetric (a normal matrix, when added with its transpose becomes symmetric)), and returns a tuple with
    both the laplacian of the matrix, and the infinity norm of the matrix.

    Example:
    J=sorted(list(gb.Matrix.ssget(2760)), key=itemgetter(0))[0]
    h=J[1]+J[1].transpose()
    L=laplacian(h)
    print(L)

    The first line is simply the way to get a matrix from the SuiteSparse collection in python using the    id of a matrix. It requires the import of "itemgetter" from operator, and the prior hand installation
    of ssgetpy.

    The second line simply makes sure that the matrix is symmetric by adding it with its transpose.

    The third line calls this laplacian function on the matrix and returns a tuple with the laplacian and the inf norm.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
*/

LAGRAPH_PUBLIC
int LAGraph_Laplacian   // compute the Laplacian matrix
(
    //Question/TODO:Should i change matricies from pointers to normal variable names?
    // outputs:
    GrB_Matrix *Lap,    // the output Laplacian matrix
    float inform,      // infinity norm of Lap
    // inputs:
    GrB_Matrix G,    // input matrix, symmetric
//  GrB_Type type,      // the type of Lap, typically GrB_FP32, ...
    char *msg
);
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

LAGRAPH_PUBLIC
int LAGraph_mypcg2
(
    //outputs
    GrB_Vector *steper,
    GrB_Index *k_result,
    // inputs:
    GrB_Matrix L,    // input matrix, symmetric, result from Laplacian
    GrB_Vector u, //vector u will be passed into another function to create Householder reflection
    float malpha,  //This float 
    GrB_Matrix invdiag,
    GrB_Vector b,
    float tol,
    float maxit, 
    //error msging
    char *msg
) ;


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

    The first line is simply the way to get a matrix from the SuiteSparse collection in python using the    id of a matrix. It requires the import of "itemgetter" from operator, and the prior hand installation

    hdip2 from the example contains the vector as the first element in the tuple, the lambda value as the second element
    with the number of iterations as a tuple with a tuple in the third element.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
*/

LAGRAPH_PUBLIC
int LAGraph_Hdip_Fiedler   // compute the Hdip_Fiedler
(
    // outputs:
    GrB_Vector *iters, //This is a vector with number of inner and outer iterations listed in that order.
    float *lamb,    // Lambda of hdip_fiedler
    GrB_Vector *x, // the hdip fielder result vector
    // inputs:
    GrB_Matrix L,    // input matrix, symmetric, result from Laplacian
    float InfNorm,
    GrB_Vector kmax,
    float emax,
    float tol,
//  GrB_Type type,      // the type of Lap, typically GrB_FP32, ...
    char *msg
);
//------------------------------------------------------------------------------
// a simple example of an algorithm
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC
GrB_Info LAGraph_HelloWorld // a simple algorithm, just for illustration
(
    // output
    GrB_Matrix *Yhandle,    // Y, created on output
    // input: not modified
    LAGraph_Graph G,
    char *msg
) ;

#endif
