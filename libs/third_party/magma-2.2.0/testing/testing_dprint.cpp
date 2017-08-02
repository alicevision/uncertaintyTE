/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016
  
       @generated from testing/testing_zprint.cpp, normal z -> d, Thu Jul 27 17:27:55 2017
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"

#if defined(__unix__) || defined(__APPLE__)
#define REDIRECT
#include <unistd.h>
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dprint
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    double *hA;
    magmaDouble_ptr dA;
    //magma_int_t ione     = 1;
    //magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t M, N, lda, ldda;  //size
    int status = 0;
    
    magma_opts opts;
    opts.parse_opts( argc, argv );

    #ifdef REDIRECT
        // dup/dup2 aren't available on Windows to restore stdout
        // save stdout and redirect to file
        const char* fname = "testing_dprint.out";
        printf( "redirecting output to %s\n", fname );
        fflush( stdout );
        int stdout_save = dup( fileno(stdout) );
        FILE* f = freopen( fname, "w", stdout );
        TESTING_CHECK( f == NULL );
    #endif

    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M     = opts.msize[itest];
            N     = opts.nsize[itest];
            lda   = M;
            ldda  = magma_roundup( M, opts.align );  // multiple of 32 by default
            //size  = lda*N;

            /* Allocate host memory for the matrix */
            TESTING_CHECK( magma_dmalloc_cpu( &hA, lda *N ));
            TESTING_CHECK( magma_dmalloc( &dA, ldda*N ));
        
            //lapackf77_dlarnv( &ione, ISEED, &size, hA );
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < M; ++i ) {
                    hA[i + j*lda] = MAGMA_D_MAKE( i + j*0.01, 0. );
                }
            }
            magma_dsetmatrix( M, N, hA, lda, dA, ldda, opts.queue );
            
            printf( "A=" );
            magma_dprint( M, N, hA, lda );
            
            printf( "dA=" );
            magma_dprint_gpu( M, N, dA, ldda, opts.queue );
            
            magma_free_cpu( hA );
            magma_free( dA );
        }
    }

    #ifdef REDIRECT
        // restore stdout
        fflush( stdout );
        dup2( stdout_save, fileno(stdout) );
        close( stdout_save );

        // compare output file to reference
        printf( "diff testing_dprint.ref testing_dprint.out\n" );
        fflush( stdout );

        int err = system( "diff testing_dprint.ref testing_dprint.out" );
        bool okay = (err == 0);
        status += ! okay;
        printf( "diff %s\n", (okay ? "ok" : "failed") );
    #endif

    opts.cleanup();
    TESTING_CHECK( magma_finalize() );
    return status;
}
