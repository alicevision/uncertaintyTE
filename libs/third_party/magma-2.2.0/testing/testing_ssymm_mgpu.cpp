/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from testing/testing_zhemm_mgpu.cpp, normal z -> s, Thu Jul 27 17:27:34 2017
       
       @author Mark Gates
       @author Azzam Haidar
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <cuda_runtime.h>

// includes, project
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_operators.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing magma_ssymm_mgpu
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    const float c_neg_one = MAGMA_S_NEG_ONE;
    const float alpha     = MAGMA_S_MAKE( 3.456, 5.678 );
    const float beta      = MAGMA_S_MAKE( 1.234, 2.456 );
    const magma_int_t ione = 1;
    
    real_Double_t    gflops, gpu_perf=0., cpu_perf=0., gpu_time=0., cpu_time=0.;
    real_Double_t    gpu_perf2=0., gpu_time2=0.;
    float           error, work[1];
    float *hA, *hB, *hC, *hR;
    magmaFloat_ptr dA[MagmaMaxGPUs], dB[MagmaMaxGPUs], dC[MagmaMaxGPUs], dwork[MagmaMaxGPUs];
    magmaFloat_ptr dA2;
    magma_int_t i, j, dev, M, N, K, Moff, Noff, Koff, size, lda, ldb, ldc, ldda, lddb, lddc, nb;
    magma_int_t iseed[4] = {0,0,0,1};
    int status = 0;
    
    magma_opts opts;
    opts.parse_opts( argc, argv );
    opts.ngpu = abs( opts.ngpu );  // always uses multi-GPU code
    
    // See testing_sgemm about tolerance.
    float eps = lapackf77_slamch("E");
    float tol = 3*eps;

    // default values
    nb = (opts.nb > 0 ? opts.nb : 64);
    
    magma_int_t nodes[MagmaMaxGPUs][MagmaMaxGPUs+2];
    magma_int_t nnode = 0;
    magma_buildconnection_mgpu( nodes, &nnode, opts.ngpu );
    
    printf("%% GPU communication pattern: nnode %lld\n", (long long) nnode);
    for (i = 0; i < nnode; ++i) {
        magma_int_t my_ndev = nodes[i][MagmaMaxGPUs];
        printf("%% node %lld has %lld GPUs:", (long long) i, (long long) my_ndev );
        for (j=0; j < my_ndev; ++j) {
            printf(" %lld", (long long) nodes[i][j] );
            if (j < my_ndev-1) {
                printf(",");
            }
        }
        printf("\n");
    }

    // number of queues per GPU. Requires ngpu.
    magma_int_t nqueue  = opts.ngpu;
    // number of events per GPU. Requires ngpu.
    magma_int_t nevents = opts.ngpu;
    magma_queue_t queues[MagmaMaxGPUs][20], queues0[MagmaMaxGPUs];
    magma_event_t events[MagmaMaxGPUs][MagmaMaxGPUs*MagmaMaxGPUs + 10];
    for( dev = 0; dev < opts.ngpu; ++dev ) {
        magma_setdevice( dev );
        for( i = 0; i < nqueue; ++i ) {
            magma_queue_create( dev, &queues[dev][i] );
        }
        queues0[dev] = queues[dev][0];
        for( i = 0; i < nevents; ++i ) {
            magma_event_create( &events[dev][i] );
        }
    }

    printf("%% nb %lld, ngpu %lld, version %lld\n", (long long) nb, (long long) opts.ngpu, (long long) opts.version );
    printf("%%   M     N    nb offset  CPU Gflop/s (sec)   GPU Gflop/s (sec)   CUBLAS hemm (sec)   ||R|| / ||A||*||B||\n");
    printf("%%========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M    = opts.msize[itest];
            N    = opts.nsize[itest];
            K    = (opts.side == MagmaLeft ? M : N);
            lda  = K;
            ldb  = M;
            ldc  = M;
            ldda = magma_roundup( lda, opts.align );  // multiple of 32 by default
            lddb = magma_roundup( ldb, opts.align );  // multiple of 32 by default
            lddc = magma_roundup( ldc, opts.align );  // multiple of 32 by default
            
            magma_int_t lwork = lddc*N + (M*N)*opts.ngpu;
            
            TESTING_CHECK( magma_smalloc_cpu( &hA, lda*K ));
            TESTING_CHECK( magma_smalloc_cpu( &hB, ldb*N ));
            TESTING_CHECK( magma_smalloc_cpu( &hC, ldc*N ));
            
            TESTING_CHECK( magma_smalloc_pinned( &hR, ldc*N ));
            
            for( dev = 0; dev < opts.ngpu; ++dev ) {
                magma_int_t Klocal = ((K / nb) / opts.ngpu + 1) * nb;  // over estimate
                magma_setdevice( dev );
                TESTING_CHECK( magma_smalloc( &dA[dev],    ldda*Klocal ));
                TESTING_CHECK( magma_smalloc( &dB[dev],    lddb*N      ));
                TESTING_CHECK( magma_smalloc( &dC[dev],    lddc*N      ));
                TESTING_CHECK( magma_smalloc( &dwork[dev], lwork       ));
            }
            
            if ( opts.check ) {
                // for running cuBLAS
                magma_setdevice( 0 );
                TESTING_CHECK( magma_smalloc( &dA2, ldda*K ));
            }
            
            for( int offset = 0; offset < N; offset += nb/2 ) {
                if (opts.side == MagmaLeft) {
                    Moff = M - offset;
                    Noff = N;
                }
                else {
                    Moff = M;
                    Noff = N - offset;
                }
                Koff = K - offset;
                gflops = FLOPS_SSYMM( MagmaLeft, Moff, Noff ) / 1e9;
                
                size = lda*K;
                lapackf77_slarnv( &ione, iseed, &size, hA );
                // set opposite triangle to NAN to ensure we don't use it.
                magma_int_t K1 = K - 1;
                if (opts.uplo == MagmaLower) 
                    lapackf77_slaset( "Upper", &K1, &K1, &MAGMA_S_NAN, &MAGMA_S_NAN, hA + lda, &lda );
                else
                    lapackf77_slaset( "Lower", &K1, &K1, &MAGMA_S_NAN, &MAGMA_S_NAN, hA + 1,   &lda );
                
                size = ldb*N;
                lapackf77_slarnv( &ione, iseed, &size, hB );
                size = ldc*N;
                lapackf77_slarnv( &ione, iseed, &size, hC );
                lapackf77_slacpy( "Full", &M, &N, hC, &ldc, hR, &ldc );
                
                // for error checks
                float Anorm = safe_lapackf77_slansy( "F", lapack_uplo_const(opts.uplo), &Koff, &hA[offset + offset*lda], &lda, work );
                float Bnorm = lapackf77_slange( "F", &Moff, &Noff, hB, &ldb, work );
                float Cnorm = lapackf77_slange( "F", &Moff, &Noff, hC, &ldc, work );
                
                /* ====================================================================
                   Performs operation using MAGMA
                   =================================================================== */
                magma_ssetmatrix_1D_col_bcyclic( opts.ngpu, K, K, nb, hA, lda, dA, ldda, queues0 );
                for( dev = 0; dev < opts.ngpu; ++dev ) {
                    magma_setdevice( dev );
                    magma_ssetmatrix( M, N, hB, ldb, dB[dev], lddb, opts.queue );
                    // send C to GPU corresponding to offset; not needed everywhere
                    magma_int_t start_dev = (offset / nb) % opts.ngpu;
                    if (dev == start_dev) {
                        magma_ssetmatrix( M, N, hC, ldc, dC[dev], lddc, opts.queue );
                    }
                }
                
                // ssymm_mgpu is synchronous, unlike async, single GPU BLAS
                gpu_time = magma_wtime();
                
                // only Left, Lower supported
                magmablas_ssymm_mgpu(
                    opts.side, opts.uplo, Moff, Noff,
                    alpha, dA, ldda, offset,
                           dB, lddb,
                    beta,  dC, lddc, dwork, lwork,
                    opts.ngpu, nb, queues, nqueue, events, nevents, nodes, nnode );
                
                gpu_time = magma_wtime() - gpu_time;
                gpu_perf = gflops / gpu_time;
                
                /* ====================================================================
                   Performs operation using CUBLAS
                   =================================================================== */
                if ( opts.check && iter == 0 ) {
                    magma_setdevice( 0 );
                    magma_ssetmatrix( K, K, hA, lda, dA2,      ldda, opts.queue );
                    magma_ssetmatrix( M, N, hB, ldb, dB[0],    lddb, opts.queue );
                    magma_ssetmatrix( M, N, hC, ldc, dwork[0], lddc, opts.queue );
                    
                    gpu_time2 = magma_sync_wtime(0);
                    magma_ssymm(
                        opts.side, opts.uplo, Moff, Noff,
                        alpha, dA2 + offset + offset*ldda, ldda,
                               dB[0],    lddb,
                        beta,  dwork[0], lddc, opts.queue );
                    gpu_time2 = magma_sync_wtime(0) - gpu_time2;
                    gpu_perf2 = gflops / gpu_time2;
                }
                
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                if ( opts.check ) {
                    
                    cpu_time = magma_wtime();
                    blasf77_ssymm( lapack_side_const(opts.side),
                                   lapack_uplo_const(opts.uplo), &Moff, &Noff,
                                   &alpha, hA + offset + offset*lda, &lda,
                                           hB, &ldb,
                                   &beta,  hC, &ldc );
                    cpu_time = magma_wtime() - cpu_time;
                    cpu_perf = gflops / cpu_time;
                    
                    for (dev=0; dev < opts.ngpu; ++dev) {
                        magma_setdevice( dev );
                        magma_sgetmatrix( M, N, dC[dev], lddc, hR, ldc, opts.queue );
                        
                        // See testing_sgemm for formula.
                        size = ldc*N;
                        blasf77_saxpy( &size, &c_neg_one, hC, &ione, hR, &ione );
                        error = lapackf77_slange( "F", &Moff, &Noff, hR, &ldc, work )
                              / (sqrt(float(Koff+2))*fabs(alpha)*Anorm*Bnorm + 2*fabs(beta)*Cnorm);
                        
                        bool okay = (error < tol);
                        status += ! okay;
                        if (dev == 0) {
                            printf( "%5lld %5lld %5lld %5lld   %7.1f (%7.4f)   %7.1f (%7.4f)   %7.1f (%7.4f)   %8.2e   %s\n",
                                    (long long) M, (long long) N, (long long) nb, (long long) offset,
                                    cpu_perf, cpu_time,
                                    gpu_perf, gpu_time,
                                    gpu_perf2, gpu_time2,
                                    error, (okay ? "ok" : "failed") );
                        }
                        else {
                            printf( "    dev %lld %74s  %8.2e   %s\n",
                                    (long long) dev, "",
                                    error, (okay ? "ok" : "failed") );
                        }
                    }
                }
                else {
                    printf( "%5lld %5lld %5lld %5lld     ---   (  ---  )   %7.1f (%7.4f)     ---   (  ---  )   ---\n",
                            (long long) M, (long long) N, (long long) nb, (long long) offset,
                            gpu_perf, gpu_time );
                }
            } // offset
            printf( "\n" );
              
            magma_free_cpu( hA );
            magma_free_cpu( hB );
            magma_free_cpu( hC );
            
            magma_free_pinned( hR );
            
            for( dev = 0; dev < opts.ngpu; ++dev ) {
                magma_setdevice( dev );
                magma_free( dA[dev]    );
                magma_free( dB[dev]    );
                magma_free( dC[dev]    );
                magma_free( dwork[dev] );
            }
            
            if ( opts.check ) {
                magma_setdevice( 0 );
                magma_free( dA2 );
            }
            fflush( stdout );
        }  // iter
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    for( dev = 0; dev < opts.ngpu; ++dev ) {
        magma_setdevice( dev );
        for( i = 0; i < nqueue; ++i ) {
            magma_queue_destroy( queues[dev][i] );
        }
        for( i = 0; i < nevents; ++i ) {
            magma_event_destroy( events[dev][i] );
        }
    }
    
    opts.cleanup();
    TESTING_CHECK( magma_finalize() );
    return status;
}
