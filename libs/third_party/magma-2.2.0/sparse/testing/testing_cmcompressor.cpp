/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from sparse/testing/testing_zmcompressor.cpp, normal z -> c, Sun Nov 20 20:20:46 2016
       @author Hartwig Anzt
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_v2.h"
#include "magmasparse.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- testing any solver
*/
int main(  int argc, char** argv )
{
    magma_int_t info = 0;
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    magma_copts zopts;
    magma_queue_t queue=NULL;
    magma_queue_create( 0, &queue );

    real_Double_t res;
    magma_c_matrix A={Magma_CSR}, AT={Magma_CSR}, A2={Magma_CSR}, 
    B={Magma_CSR}, dB={Magma_CSR};
    
    int i=1;
    real_Double_t start, end;
    TESTING_CHECK( magma_cparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while( i < argc ) {
        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_int_t laplace_size = atoi( argv[i] );
            TESTING_CHECK( magma_cm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            TESTING_CHECK( magma_c_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "\n# matrix info: %lld-by-%lld with %lld nonzeros\n\n",
                (long long) A.num_rows, (long long) A.num_cols, (long long) A.nnz );

        // scale matrix
        TESTING_CHECK( magma_cmscale( &A, zopts.scaling, queue ));

        // remove nonzeros in matrix
        start = magma_sync_wtime( queue );
        for (int j=0; j < 10; j++) {
            TESTING_CHECK( magma_cmcsrcompressor( &A, queue ));
        }
        end = magma_sync_wtime( queue );
        printf( " > MAGMA CPU: %.2e seconds.\n", (end-start)/10 );
        // transpose
        TESTING_CHECK( magma_cmtranspose( A, &AT, queue ));

        // convert, copy back and forth to check everything works
        TESTING_CHECK( magma_cmconvert( AT, &B, Magma_CSR, Magma_CSR, queue ));
        magma_cmfree(&AT, queue );
        TESTING_CHECK( magma_cmtransfer( B, &dB, Magma_CPU, Magma_DEV, queue ));
        magma_cmfree(&B, queue );

        start = magma_sync_wtime( queue );
        for (int j=0; j < 10; j++) {
            TESTING_CHECK( magma_cmcsrcompressor_gpu( &dB, queue ));
        }
        end = magma_sync_wtime( queue );
        printf( " > MAGMA GPU: %.2e seconds.\n", (end-start)/10 );


        TESTING_CHECK( magma_cmtransfer( dB, &B, Magma_DEV, Magma_CPU, queue ));
        magma_cmfree(&dB, queue );
        TESTING_CHECK( magma_cmconvert( B, &AT, Magma_CSR, Magma_CSR, queue ));
        magma_cmfree(&B, queue );

        // transpose back
        TESTING_CHECK( magma_cmtranspose( AT, &A2, queue ));
        magma_cmfree(&AT, queue );
        TESTING_CHECK( magma_cmdiff( A, A2, &res, queue ));
        printf("%% ||A-B||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("%% tester matrix compressor:  ok\n");
        else
            printf("%% tester matrix compressor:  failed\n");

        magma_cmfree(&A, queue );
        magma_cmfree(&A2, queue );

        i++;
    }
    
    magma_queue_destroy( queue );
    TESTING_CHECK( magma_finalize() );
    return info;
}
