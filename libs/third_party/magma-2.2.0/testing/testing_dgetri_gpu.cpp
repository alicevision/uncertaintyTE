/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from testing/testing_zgetri_gpu.cpp, normal z -> d, Thu Jul 27 17:28:26 2017
       @author Mark Gates
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
   -- Testing dgetri
*/
int main( int argc, char** argv )
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    // constants
    const double c_zero    = MAGMA_D_ZERO;
    const double c_one     = MAGMA_D_ONE;
    const double c_neg_one = MAGMA_D_NEG_ONE;
    
    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double *h_A, *h_Ainv, *h_R, *work, unused[1];
    magmaDouble_ptr d_A, dwork;
    magma_int_t N, n2, lda, ldda, info, lwork, ldwork;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    double tmp;
    double error, rwork[1];
    magma_int_t *ipiv, iunused[1];
    int status = 0;
    
    magma_opts opts;
    opts.parse_opts( argc, argv );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("%%   N   CPU Gflop/s (sec)   GPU Gflop/s (sec)   ||I - A*A^{-1}||_1 / (N*cond(A))\n");
    printf("%%===============================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            ldda   = magma_roundup( N, opts.align );  // multiple of 32 by default
            ldwork = N * magma_get_dgetri_nb( N );
            gflops = FLOPS_DGETRI( N ) / 1e9;
            
            // query for workspace size
            lwork = -1;
            lapackf77_dgetri( &N, unused, &lda, iunused, &tmp, &lwork, &info );
            if (info != 0) {
                printf("lapackf77_dgetri returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            lwork = magma_int_t( MAGMA_D_REAL( tmp ));
            
            TESTING_CHECK( magma_imalloc_cpu( &ipiv,   N      ));
            TESTING_CHECK( magma_dmalloc_cpu( &work,   lwork  ));
            TESTING_CHECK( magma_dmalloc_cpu( &h_A,    n2     ));
            TESTING_CHECK( magma_dmalloc_cpu( &h_Ainv, n2     ));
            TESTING_CHECK( magma_dmalloc_cpu( &h_R,    n2     ));
            
            TESTING_CHECK( magma_dmalloc( &d_A,    ldda*N ));
            TESTING_CHECK( magma_dmalloc( &dwork,  ldwork ));
            
            /* Initialize the matrix */
            lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
            
            /* Factor the matrix. Both MAGMA and LAPACK will use this factor. */
            magma_dsetmatrix( N, N, h_A, lda, d_A, ldda, opts.queue );
            magma_dgetrf_gpu( N, N, d_A, ldda, ipiv, &info );
            magma_dgetmatrix( N, N, d_A, ldda, h_Ainv, lda, opts.queue );
            if (info != 0) {
                printf("magma_dgetrf_gpu returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            // check for exact singularity
            //h_Ainv[ 10 + 10*lda ] = MAGMA_D_MAKE( 0.0, 0.0 );
            //magma_dsetmatrix( N, N, h_Ainv, lda, d_A, ldda, opts.queue );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_wtime();
            magma_dgetri_gpu( N, d_A, ldda, ipiv, dwork, ldwork, &info );
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0) {
                printf("magma_dgetri_gpu returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                lapackf77_dgetri( &N, h_Ainv, &lda, ipiv, work, &lwork, &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0) {
                    printf("lapackf77_dgetri returned error %lld: %s.\n",
                           (long long) info, magma_strerror( info ));
                }
                printf( "%5lld   %7.2f (%7.2f)   %7.2f (%7.2f)",
                        (long long) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf( "%5lld     ---   (  ---  )   %7.2f (%7.2f)",
                        (long long) N, gpu_perf, gpu_time );
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.check ) {
                magma_dgetmatrix( N, N, d_A, ldda, h_Ainv, lda, opts.queue );
                
                // compute 1-norm condition number estimate, following LAPACK's zget03
                double normA, normAinv, rcond;
                normA    = lapackf77_dlange( "1", &N, &N, h_A,    &lda, rwork );
                normAinv = lapackf77_dlange( "1", &N, &N, h_Ainv, &lda, rwork );
                if ( normA <= 0 || normAinv <= 0 ) {
                    rcond = 0;
                    error = 1 / (tol/opts.tolerance);  // == 1/eps
                }
                else {
                    rcond = (1 / normA) / normAinv;
                    // R = I
                    // R -= A*A^{-1}
                    // err = ||I - A*A^{-1}|| / ( N ||A||*||A^{-1}|| ) = ||R|| * rcond / N, using 1-norm
                    lapackf77_dlaset( "full", &N, &N, &c_zero, &c_one, h_R, &lda );
                    blasf77_dgemm( "no", "no", &N, &N, &N,
                                   &c_neg_one, h_A,    &lda,
                                               h_Ainv, &lda,
                                   &c_one,     h_R,    &lda );
                    error = lapackf77_dlange( "1", &N, &N, h_R, &lda, rwork );
                    error = error * rcond / N;
                }
                
                bool okay = (error < tol);
                status += ! okay;
                printf( "   %8.2e   %s\n",
                        error, (okay ? "ok" : "failed"));
            }
            else {
                printf( "\n" );
            }
            
            magma_free_cpu( ipiv   );
            magma_free_cpu( work   );
            magma_free_cpu( h_A    );
            magma_free_cpu( h_Ainv );
            magma_free_cpu( h_R    );
            
            magma_free( d_A    );
            magma_free( dwork  );
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
