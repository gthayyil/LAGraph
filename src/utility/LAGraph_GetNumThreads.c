//------------------------------------------------------------------------------
// LAGraph_GetNumThreads: get the # of threads to use
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_get_nthreads: get # of threads that will be used by LAGraph.

#include "LG_internal.h"

int LAGraph_GetNumThreads
(
    int *nthreads,      // # of threads to use
    char *msg
)
{

   //---------------------------------------------------------------------------
   // check inputs
   //---------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (nthreads != NULL, GrB_NULL_POINTER) ;

   //---------------------------------------------------------------------------
   // get number of threads
   //---------------------------------------------------------------------------

    #if LG_SUITESPARSE
    {
        // SuiteSparse:GraphBLAS: get # of threads from global setting
        GrB_TRY (GxB_get (GxB_NTHREADS, nthreads)) ;
    }
    #elif defined ( _OPENMP )
    {
        // get # of threads from OpenMP global setting
        (*nthreads) = omp_get_max_threads ( ) ;
    }
    #else
    {
        // single-threaded if not using SuiteSparse:GraphBLAS or OpenMP
        (*nthreads) = 1 ;
    }
    #endif
    return (GrB_SUCCESS) ;
}

