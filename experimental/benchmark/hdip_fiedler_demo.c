//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/hdip_fielder_demo.c: a simple demo
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A Davis, Texas A&M University

//------------------------------------------------------------------------------

// This main program is a simple driver for testing and benchmarking the
// LAGraph_HelloWorld "algorithm", in experimental/algorithm.  To use it,
// compile LAGraph while in the build folder with these commands:
//
//      cd LAGraph/build
//      cmake ..
//      make -j8
//
// Then run this demo with an input matrix.  For example:
//
//      ./experimental/benchmark/hellworld_demo ../data/west0067.mtx
//      ./experimental/benchmark/hellworld_demo < ../data/west0067.mtx
//      ./experimental/benchmark/hellworld_demo ../data/karate.mtx
//
// If you create your own algorithm and want to mimic this main program, call
// it write in experimental/benchmark/whatever_demo.c (with "_demo.c" as the
// end of the filename), and the cmake will find it and compile it.

// This main program makes use of supporting utilities in
// src/benchmark/LAGraph_demo.h and src/utility/LG_internal.h.
// See helloworld2_demo.c for a main program that just uses the
// user-callable methods in LAGraph.h and LAGraphX.h.

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"

// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&Y) ;                             \
    LAGraph_Delete (&G, msg) ;                  \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix Y = NULL ;
    GrB_Matrix A = NULL;

    // start GraphBLAS and LAGraph
    bool burble = false ;               // set true for diagnostic outputs
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph: this method is defined in LAGraph_demo.h
    //--------------------------------------------------------------------------

    // readproblem can read in a file in Matrix Market format, or in a binary
    // format created by binwrite (see LAGraph_demo.h, or the main program,
    // mtx2bin_demo).

    double t = LAGraph_WallClockTime ( ) ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LG_TRY (readproblem (
        &G,         // the graph that is read from stdin or a file
        NULL,       // source nodes (none, if NULL)
        true,      // make the graph undirected, if true
        true,      // remove self-edges, if true
        false,      // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LG_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // try the LAGraph_HelloWorld "algorithm"
    //--------------------------------------------------------------------------
    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;
    LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;
    GRB_TRY (GrB_Matrix_new (&A, GrB_FP32, n, n)) ;
    GRB_TRY (GrB_assign (A, G->A, NULL, (double) 1,
        GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;
    GrB_free (&(G->A)) ;
    G->A = A ;

    GRB_TRY (GrB_Matrix_new (&Y, GrB_FP32, n, n)) ; 
    
    //Variables needed to test Laplacian
    float inform;
    LG_TRY (LAGraph_Laplacian(&Y,inform, A, msg)) ;
    
    //Additional variables and modifications needed to test MYPCG2
    GrB_Vector steper = NULL;
    GrB_Vector u = NULL; // a vector of size nrowsLap, filled with 1.
    //set u[0] = 1+sqrt(nrowsLap)
    GrB_Matrix indiag = NULL;
    GrB_Vector x = NULL;
    GrB_Index k;
    GrB_Index nrows;
    float nrowsLap; //number of rows of laplacian matrix
    float alpha;
    GRB_TRY (GrB_Matrix_nrows(&nrows,Y ));
    
    nrowsLap = nrows;
    GRB_TRY (GrB_Vector_new (&u, GrB_FP32, nrowsLap));
    GRB_TRY (GrB_apply (u, NULL, NULL, GxB_ONE_FP32, u, NULL)) ;
    GRB_TRY (GrB_Vector_setElement_FP32(u, 1+sqrt(nrowsLap), 0));
    
    alpha = nrowsLap + sqrt(nrowsLap);
    
    //printf("works?4");
    GRB_TRY (GrB_Matrix_new (&indiag, GrB_FP32, nrows,nrows));
    GRB_TRY (GrB_select (indiag, NULL, NULL, GrB_DIAG, Y,0, NULL));
    //printf("works?4");
    GRB_TRY (GrB_apply (indiag, NULL, NULL, GrB_MINV_FP32, indiag, NULL));
    
    //printf("works?3");
    GRB_TRY (GrB_Vector_new (&x, GrB_FP32, nrowsLap));
    GRB_TRY (GrB_apply (x, NULL, NULL, GxB_ONE_FP32, x, NULL));
    GRB_TRY (GrB_Vector_setElement_FP32(x, 0, 0));
    
    //printf("works?2");
    //t = LAGraph_WallClockTime( );
    printf("works?");
    LG_TRY (LAGraph_mypcg2(steper, k, Y, u, alpha, indiag, x, .000001, 50, msg));
    printf("aftermypcg2");
    //t = LAGraph_WallClockTime( ) - t;
    printf ("Time for LAGraph_HelloWorld: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // check the results (make sure Y is a copy of G->A)
    //--------------------------------------------------------------------------

    bool isequal ;
    t = LAGraph_WallClockTime ( ) ;
    //LG_TRY (LAGraph_Matrix_IsEqual (&isequal, Y, G->A, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check results:       %g sec\n", t) ;
    if (isequal)
    {
        printf ("Test passed.\n") ;
    }
    else
    {
        printf ("Test failure!\n") ;
    }

    //--------------------------------------------------------------------------
    // print the results (Y is just a copy of G->A)
    //--------------------------------------------------------------------------

    printf ("\n===============================The result matrix Y:\n") ;
    LG_TRY (LAGraph_Matrix_Print (Y, LAGraph_SHORT, stdout, msg)) ;
    LG_TRY (LAGraph_Vector_Print (steper, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------
         
                          
                        
    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

