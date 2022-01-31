//------------------------------------------------------------------------------
// LAGraph_SSaveSet: save a set of matrices to a *.lagraph file
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// Contributed by Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAGraph_SSaveSet saves a set of matrices to a *.lagraph file.
// The file is created, written to with the JSON header and the serialized
// matrices, and then closed.  If using SuiteSparse:GraphBLAS, the highest
// level of compression is used (LZ4HC:9).

// Use LAGraph_SSLoadSet to load the matrices back in from the file.

// This method will not work without SuiteSparse:GraphBLAS, because the C API
// has no GrB* method for querying the GrB_Type (or its name as a string) of a
// matrix.

//------------------------------------------------------------------------------

#define LAGraph_FREE_WORK                           \
{                                                   \
    fclose (f) ;                                    \
    f = NULL ;                                      \
    GrB_free (&desc) ;                              \
    LAGraph_SFreeContents (&Contents, nmatrices) ;  \
}

#define LAGraph_FREE_ALL                            \
{                                                   \
    LAGraph_FREE_WORK ;                             \
}

#include "LG_internal.h"
#include "LAGraphX.h"

//------------------------------------------------------------------------------
// LAGraph_SSaveSet
//------------------------------------------------------------------------------

int LAGraph_SSaveSet            // save a set of matrices from a *.lagraph file
(
    // inputs:
    char *filename,             // name of file to write to
    GrB_Matrix *Set,            // array of GrB_Matrix of size nmatrices
    GrB_Index nmatrices,        // # of matrices to write to *.lagraph file
    char *collection,           // name of this collection of matrices
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    FILE *f = NULL ;

    LAGraph_Contents *Contents = NULL ;
    GrB_Descriptor desc = NULL ;

    LG_ASSERT (filename != NULL && Set != NULL && collection != NULL,
        GrB_NULL_POINTER) ;

    #if LG_SUITESPARSE
    GrB_TRY (GrB_Descriptor_new (&desc)) ;
    GrB_TRY (GxB_set (desc, GxB_COMPRESSION, GxB_COMPRESSION_LZ4HC + 9)) ;
    #endif

    f = fopen (filename, "w") ;
    LG_ASSERT_MSG (f != NULL, -1001, "unable to create output file") ;

    //--------------------------------------------------------------------------
    // serialize all the matrices
    //--------------------------------------------------------------------------

    // allocate an Contents array of size nmatrices to hold the contents
    Contents = LAGraph_Calloc (nmatrices, sizeof (LAGraph_Contents)) ;
    LG_ASSERT (Contents != NULL, GrB_OUT_OF_MEMORY) ;

    for (GrB_Index i = 0 ; i < nmatrices ; i++)
    {
        #if LG_SUITESPARSE
        {
            GrB_TRY (GxB_Matrix_serialize (&(Contents [i].blob),
                (GrB_Index *)&(Contents [i].blob_size), Set [i], desc)) ;
        }
        #else
        {
            GrB_Index estimate ;
            GrB_TRY (GrB_Matrix_serializeSize (&estimate, Set [i])) ;
            Contents [i].blob_size = estimate ;
            Contents [i].blob = LAGraph_Malloc (estimate, sizeof (uint8_t)) ;
            LG_ASSERT (Contents [i].blob != NULL, GrB_OUT_OF_MEMORY) ;
            GrB_TRY (GrB_Matrix_serialize (Contents [i].blob,
                (GrB_Index *) &(Contents [i].blob_size), Set [i])) ;
            bool ok ;
            Contents [i].blob = LAGraph_Realloc (&(Contents [i].blob_size),
                estimate, sizeof (uint8_t), Contents [i].blob, &ok) ;
            LG_ASSERT (ok, GrB_OUT_OF_MEMORY) ;
        }
        #endif
    }

    //--------------------------------------------------------------------------
    // write the header
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_SWrite_HeaderStart (f, collection, msg)) ;
    for (GrB_Index i = 0 ; i < nmatrices ; i++)
    {
        #if LG_SUITESPARSE
        char typename [GxB_MAX_NAME_LEN] ;
        GrB_TRY (GxB_Matrix_type_name (typename, Set [i])) ;
        #else
        // This will fail:  the C API has no method for querying the type of a
        // matrix, so SuiteSparse is required for this LAGraph_SSaveSet to
        // work. The C API urgently needs GrB_Matrix_type_name added to it.
        char *typename = NULL ;
        #endif
        char matrix_name [256] ;
        snprintf (matrix_name, 256, "A_%" PRIu64, i) ;
        LG_TRY (LAGraph_SWrite_HeaderItem (f, LAGraph_matrix_kind,
            matrix_name, typename, 0, Contents [i].blob_size, msg)) ;
    }
    LG_TRY (LAGraph_SWrite_HeaderEnd (f, msg)) ;

    //--------------------------------------------------------------------------
    // write all the blobs
    //--------------------------------------------------------------------------

    for (GrB_Index i = 0 ; i < nmatrices ; i++)
    {
        LG_TRY (LAGraph_SWrite_Item (f, Contents [i].blob,
            Contents [i].blob_size, msg)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

