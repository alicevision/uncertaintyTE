/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @author Mark Gates
       @generated from testing/testing_zlag2c.cpp, mixed zc -> ds, Thu Jul 27 17:27:44 2017
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// includes, project
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dlag2s and slag2d
*/
int main( int argc, char** argv )
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();
    
    real_Double_t   gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double error, work[1];
    float serror, swork[1];
    double c_neg_one = MAGMA_D_NEG_ONE;
    float  s_neg_one = MAGMA_S_NEG_ONE;
    magma_int_t ione = 1;
    magma_int_t m, n, lda, ldda, size, info;
    magma_int_t ISEED[4] = {0,0,0,1};
    int status = 0;
    float   *SA, *SR;
    double   *A,  *R;
    magmaFloat_ptr dSA;
    magmaDouble_ptr dA;
    
    magma_opts opts;
    opts.parse_opts( argc, argv );
    
    printf("%% func     M     N     CPU GB/s (ms)       GPU GB/s (ms)     ||R||_F\n");
    printf("%%====================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            m = opts.msize[itest];
            n = opts.nsize[itest];
            lda  = m;
            ldda = magma_roundup( m, opts.align );  // multiple of 32 by default
            // m*n double-real loads and m*n single-real stores (and vice-versa for slag2d)
            gbytes = (real_Double_t) m*n * (sizeof(double) + sizeof(float)) / 1e9;
            size = ldda*n;  // ldda >= lda
            
            TESTING_CHECK( magma_smalloc_cpu( &SA, size ));
            TESTING_CHECK( magma_dmalloc_cpu( &A, size ));
            TESTING_CHECK( magma_smalloc_cpu( &SR, size ));
            TESTING_CHECK( magma_dmalloc_cpu( &R, size ));
            
            TESTING_CHECK( magma_smalloc( &dSA, size ));
            TESTING_CHECK( magma_dmalloc( &dA, size ));
            
            lapackf77_dlarnv( &ione, ISEED, &size,  A );
            lapackf77_slarnv( &ione, ISEED, &size, SA );
            
            magma_dsetmatrix( m, n, A,  lda, dA,  ldda, opts.queue );
            magma_ssetmatrix( m, n, SA, lda, dSA, ldda, opts.queue );
            
            /* =====================================================================
               Performs operation using LAPACK dlag2s
               =================================================================== */
            cpu_time = magma_wtime();
            lapackf77_dlag2s( &m, &n, A, &lda, SA, &lda, &info );
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (info != 0) {
                printf("lapackf77_dlag2s returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            /* ====================================================================
               Performs operation using MAGMA dlag2s
               =================================================================== */
            gpu_time = magma_sync_wtime( opts.queue );
            magmablas_dlag2s( m, n, dA, ldda, dSA, ldda, opts.queue, &info );
            gpu_time = magma_sync_wtime( opts.queue ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (info != 0) {
                printf("magmablas_dlag2s returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            magma_sgetmatrix( m, n, dSA, ldda, SR, lda, opts.queue );
            
            /* =====================================================================
               compute error |SA_magma - SA_lapack|
               should be zero if both are IEEE compliant
               =================================================================== */
            blasf77_saxpy( &size, &s_neg_one, SA, &ione, SR, &ione );
            serror = lapackf77_slange( "Fro", &m, &n, SR, &lda, swork );
            
            printf( "dlag2s %5lld %5lld   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                    (long long) m, (long long) n,
                    cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                    serror, (serror == 0 ? "ok" : "failed") );
            status += ! (serror == 0);
            
            /* =====================================================================
               Reset matrices
               =================================================================== */
            lapackf77_dlarnv( &ione, ISEED, &size,  A );
            lapackf77_slarnv( &ione, ISEED, &size, SA );
            
            magma_dsetmatrix( m, n, A,  lda, dA,  ldda, opts.queue );
            magma_ssetmatrix( m, n, SA, lda, dSA, ldda, opts.queue );
            
            /* =====================================================================
               Performs operation using LAPACK slag2d
               =================================================================== */
            cpu_time = magma_wtime();
            lapackf77_slag2d( &m, &n, SA, &lda, A, &lda, &info );
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (info != 0) {
                printf("lapackf77_slag2d returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            /* ====================================================================
               Performs operation using MAGMA slag2d
               =================================================================== */
            magma_ssetmatrix( m, n, SA, lda, dSA, ldda, opts.queue );
            
            gpu_time = magma_sync_wtime( opts.queue );
            magmablas_slag2d( m, n, dSA, ldda, dA, ldda, opts.queue, &info );
            gpu_time = magma_sync_wtime( opts.queue ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (info != 0) {
                printf("magmablas_slag2d returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            magma_dgetmatrix( m, n, dA, ldda, R, lda, opts.queue );
            
            /* =====================================================================
               compute error |A_magma - A_lapack|
               should be zero if both are IEEE compliant
               =================================================================== */
            blasf77_daxpy( &size, &c_neg_one, A, &ione, R, &ione );
            error = lapackf77_dlange( "Fro", &m, &n, R, &lda, work );
            
            printf( "slag2d %5lld %5lld   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                    (long long) m, (long long) n,
                    cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                    error, (error == 0 ? "ok" : "failed") );
            status += ! (error == 0);
            
            magma_free_cpu( SA );
            magma_free_cpu( A );
            magma_free_cpu( SR );
            magma_free_cpu( R );
            
            magma_free( dSA );
            magma_free( dA );
            printf( "\n" );
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
