/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from testing/testing_zgels3_gpu.cpp, normal z -> s, Thu Jul 27 17:28:37 2017

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgels
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           gpu_error, cpu_error, error, Anorm, work[1];
    float  c_one     = MAGMA_S_ONE;
    float  c_neg_one = MAGMA_S_NEG_ONE;
    float *h_A, *h_A2, *h_B, *h_X, *h_R, *tau, *h_work, tmp[1], unused[1];
    magmaFloat_ptr d_A, d_B;
    magma_int_t M, N, size, nrhs, lda, ldb, ldda, lddb, min_mn, max_mn, nb, info;
    magma_int_t lworkgpu, lhwork, lhwork2;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    magma_opts opts;
    opts.parse_opts( argc, argv );
 
    int status = 0;
    float tol = opts.tolerance * lapackf77_slamch("E");

    nrhs = opts.nrhs;
    
    printf("%%                                                           ||b-Ax|| / (N||A||)   ||dx-x||/(N||A||)\n");
    printf("%%   M     N  NRHS   CPU Gflop/s (sec)   GPU Gflop/s (sec)   CPU        GPU                         \n");
    printf("%%==================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            if ( M < N ) {
                printf( "%5lld %5lld %5lld   skipping because M < N is not yet supported.\n", (long long) M, (long long) N, (long long) nrhs );
                continue;
            }
            min_mn = min(M, N);
            max_mn = max(M, N);
            lda    = M;
            ldb    = max_mn;
            size   = lda*N;
            ldda   = magma_roundup( M, opts.align );  // multiple of 32 by default
            lddb   = magma_roundup( max_mn, opts.align );  // multiple of 32 by default
            nb     = magma_get_sgeqrf_nb( M, N );
            gflops = (FLOPS_SGEQRF( M, N ) + FLOPS_SGEQRS( M, N, nrhs )) / 1e9;
            
            lworkgpu = (M - N + nb)*(nrhs + nb) + nrhs*nb;
            
            // query for workspace size
            lhwork = -1;
            lapackf77_sgeqrf( &M, &N, unused, &M, unused, tmp, &lhwork, &info );
            lhwork2 = (magma_int_t) MAGMA_S_REAL( tmp[0] );
            
            lhwork = -1;
            lapackf77_sormqr( MagmaLeftStr, MagmaTransStr,
                              &M, &nrhs, &min_mn,
                              unused, &lda, unused,
                              unused, &ldb,
                              tmp, &lhwork, &info);
            lhwork = (magma_int_t) MAGMA_S_REAL( tmp[0] );
            lhwork = max( max( lhwork, lhwork2 ), lworkgpu );
            
            TESTING_CHECK( magma_smalloc_cpu( &tau,    min_mn    ));
            TESTING_CHECK( magma_smalloc_cpu( &h_A,    lda*N     ));
            TESTING_CHECK( magma_smalloc_cpu( &h_A2,   lda*N     ));
            TESTING_CHECK( magma_smalloc_cpu( &h_B,    ldb*nrhs  ));
            TESTING_CHECK( magma_smalloc_cpu( &h_X,    ldb*nrhs  ));
            TESTING_CHECK( magma_smalloc_cpu( &h_R,    ldb*nrhs  ));
            TESTING_CHECK( magma_smalloc_cpu( &h_work, lhwork    ));
            
            TESTING_CHECK( magma_smalloc( &d_A,    ldda*N    ));
            TESTING_CHECK( magma_smalloc( &d_B,    lddb*nrhs ));
            
            /* Initialize the matrices */
            lapackf77_slarnv( &ione, ISEED, &size, h_A );
            lapackf77_slacpy( MagmaFullStr, &M, &N, h_A, &lda, h_A2, &lda );
            
            // make random RHS
            size = M*nrhs;
            lapackf77_slarnv( &ione, ISEED, &size, h_B );
            lapackf77_slacpy( MagmaFullStr, &M, &nrhs, h_B, &ldb, h_R, &ldb );
            
            // make consistent RHS
            //size = N*nrhs;
            //lapackf77_slarnv( &ione, ISEED, &size, h_X );
            //blasf77_sgemm( MagmaNoTransStr, MagmaNoTransStr, &M, &nrhs, &N,
            //               &c_one,  h_A, &lda,
            //                        h_X, &ldb,
            //               &c_zero, h_B, &ldb );
            //lapackf77_slacpy( MagmaFullStr, &M, &nrhs, h_B, &ldb, h_R, &ldb );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_ssetmatrix( M, N,    h_A, lda, d_A, ldda, opts.queue );
            magma_ssetmatrix( M, nrhs, h_B, ldb, d_B, lddb, opts.queue );
            
            gpu_time = magma_wtime();
            magma_sgels3_gpu( MagmaNoTrans, M, N, nrhs, d_A, ldda,
                              d_B, lddb, h_work, lworkgpu, &info);
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0) {
                printf("magma_sgels3_gpu returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            // Get the solution in h_X
            magma_sgetmatrix( N, nrhs, d_B, lddb, h_X, ldb, opts.queue );
            
            // compute the residual
            blasf77_sgemm( MagmaNoTransStr, MagmaNoTransStr, &M, &nrhs, &N,
                           &c_neg_one, h_A, &lda,
                                       h_X, &ldb,
                           &c_one,     h_R, &ldb);
            Anorm = lapackf77_slange("f", &M, &N, h_A, &lda, work);
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            lapackf77_slacpy( MagmaFullStr, &M, &nrhs, h_B, &ldb, h_X, &ldb );
            
            cpu_time = magma_wtime();
            lapackf77_sgels( MagmaNoTransStr, &M, &N, &nrhs,
                             h_A, &lda, h_X, &ldb, h_work, &lhwork, &info);
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0) {
                printf("lapackf77_sgels returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            blasf77_sgemm( MagmaNoTransStr, MagmaNoTransStr, &M, &nrhs, &N,
                           &c_neg_one, h_A2, &lda,
                                       h_X,  &ldb,
                           &c_one,     h_B,  &ldb);
            
            cpu_error = lapackf77_slange("f", &M, &nrhs, h_B, &ldb, work) / (min_mn*Anorm);
            gpu_error = lapackf77_slange("f", &M, &nrhs, h_R, &ldb, work) / (min_mn*Anorm);
            
            // error relative to LAPACK
            size = M*nrhs;
            blasf77_saxpy( &size, &c_neg_one, h_B, &ione, h_R, &ione );
            error = lapackf77_slange("f", &M, &nrhs, h_R, &ldb, work) / (min_mn*Anorm);
            
            printf("%5lld %5lld %5lld   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %8.2e",
                   (long long) M, (long long) N, (long long) nrhs,
                   cpu_perf, cpu_time, gpu_perf, gpu_time, cpu_error, gpu_error, error );
            
            if ( M == N ) {
                printf( "   %s\n", (gpu_error < tol && error < tol ? "ok" : "failed"));
                status += ! (gpu_error < tol && error < tol);
            }
            else {
                printf( "   %s\n", (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }

            magma_free_cpu( tau    );
            magma_free_cpu( h_A    );
            magma_free_cpu( h_A2   );
            magma_free_cpu( h_B    );
            magma_free_cpu( h_X    );
            magma_free_cpu( h_R    );
            magma_free_cpu( h_work );
            
            magma_free( d_A    );
            magma_free( d_B    );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    opts.cleanup();
    TESTING_CHECK( magma_finalize() );
    return status;
}
