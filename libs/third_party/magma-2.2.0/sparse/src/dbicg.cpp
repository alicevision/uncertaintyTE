/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from sparse/src/zbicg.cpp, normal z -> d, Sun Nov 20 20:20:44 2016
       @author Hartwig Anzt

*/
#include "magmasparse_internal.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a general matrix.
    This is a GPU implementation of the Biconjugate Gradient method.

    Arguments
    ---------

    @param[in]
    A           magma_d_matrix
                input matrix A

    @param[in]
    b           magma_d_matrix
                RHS b

    @param[in,out]
    x           magma_d_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_d_solver_par*
                solver parameters
    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_dgesv
    ********************************************************************/

extern "C" magma_int_t
magma_dbicg(
    magma_d_matrix A, magma_d_matrix b, magma_d_matrix *x,
    magma_d_solver_par *solver_par,
    magma_queue_t queue )
{
    magma_int_t info = MAGMA_NOTCONVERGED;
    
    // prepare solver feedback
    solver_par->solver = Magma_BICG;
    solver_par->numiter = 0;
    solver_par->spmv_count = 0;

    // some useful variables
    double c_zero = MAGMA_D_ZERO;
    double c_one  = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;
    
    magma_int_t dofs = A.num_rows * b.num_cols;

    // workspace
    magma_d_matrix r={Magma_CSR}, rt={Magma_CSR}, p={Magma_CSR}, pt={Magma_CSR}, 
                z={Magma_CSR}, zt={Magma_CSR}, q={Magma_CSR}, y={Magma_CSR}, 
                yt={Magma_CSR},  qt={Magma_CSR};
                
    // need to transpose the matrix
    magma_d_matrix AT={Magma_CSR}, Ah1={Magma_CSR}, Ah2={Magma_CSR};
    
    CHECK( magma_dvinit( &r, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &rt,Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &p, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &pt,Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &q, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &qt,Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &y, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &yt,Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &z, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &zt,Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));

    
    // solver variables
    double alpha, rho, beta, rho_new, ptq;
    double res, nomb, nom0, r0;

        // transpose the matrix
    magma_dmtransfer( A, &Ah1, Magma_DEV, Magma_CPU, queue );
    magma_dmconvert( Ah1, &Ah2, A.storage_type, Magma_CSR, queue );
    magma_dmfree(&Ah1, queue );
    magma_dmtransposeconjugate( Ah2, &Ah1, queue );
    magma_dmfree(&Ah2, queue );
    Ah2.blocksize = A.blocksize;
    Ah2.alignment = A.alignment;
    magma_dmconvert( Ah1, &Ah2, Magma_CSR, A.storage_type, queue );
    magma_dmfree(&Ah1, queue );
    magma_dmtransfer( Ah2, &AT, Magma_CPU, Magma_DEV, queue );
    magma_dmfree(&Ah2, queue );
    
    // solver setup
    CHECK(  magma_dresidualvec( A, b, *x, &r, &nom0, queue));
    res = nom0;
    solver_par->init_res = nom0;
    magma_dcopy( dofs, r.dval, 1, rt.dval, 1, queue );                  // rr = r
    rho_new = magma_ddot( dofs, rt.dval, 1, r.dval, 1, queue );             // rho=<rr,r>
    rho = alpha = MAGMA_D_MAKE( 1.0, 0. );

    nomb = magma_dnrm2( dofs, b.dval, 1, queue );
    if ( nomb == 0.0 ){
        nomb=1.0;
    }       
    if ( (r0 = nomb * solver_par->rtol) < ATOLERANCE ){
        r0 = ATOLERANCE;
    }
    
    solver_par->final_res = solver_par->init_res;
    solver_par->iter_res = solver_par->init_res;
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = nom0;
        solver_par->timing[0] = 0.0;
    }
    if ( nom0 < r0 ) {
        info = MAGMA_SUCCESS;
        goto cleanup;
    }

    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_sync_wtime( queue );


    solver_par->numiter = 0;
    solver_par->spmv_count = 0;
    // start iteration
    do
    {
        solver_par->numiter++;

        magma_dcopy( dofs, r.dval, 1 , y.dval, 1, queue );             // y=r
        magma_dcopy( dofs, y.dval, 1 , z.dval, 1, queue );             // z=y
        magma_dcopy( dofs, rt.dval, 1 , yt.dval, 1, queue );           // yt=rt
        magma_dcopy( dofs, yt.dval, 1 , zt.dval, 1, queue );           // zt=yt
        
        rho= rho_new;
        rho_new = magma_ddot( dofs, rt.dval, 1, z.dval, 1, queue );  // rho=<rt,z>
        if( magma_d_isnan_inf( rho_new ) ){
            info = MAGMA_DIVERGENCE;
            break;
        }
        
        if( solver_par->numiter==1 ){
            magma_dcopy( dofs, z.dval, 1 , p.dval, 1, queue );           // yt=rt
            magma_dcopy( dofs, zt.dval, 1 , pt.dval, 1, queue );           // zt=yt
        } else {
            beta = rho_new/rho;
            magma_dscal( dofs, beta, p.dval, 1, queue );                 // p = beta*p
            magma_daxpy( dofs, c_one , z.dval, 1 , p.dval, 1, queue );   // p = z+beta*p
            magma_dscal( dofs, MAGMA_D_CONJ(beta), pt.dval, 1, queue );   // pt = beta*pt
            magma_daxpy( dofs, c_one , zt.dval, 1 , pt.dval, 1, queue );  // pt = zt+beta*pt
        }
        CHECK( magma_d_spmv( c_one, A, p, c_zero, q, queue ));      // v = Ap
        CHECK( magma_d_spmv( c_one, AT, pt, c_zero, qt, queue ));   // v = Ap
        solver_par->spmv_count++;
        solver_par->spmv_count++;
        ptq = magma_ddot( dofs, pt.dval, 1, q.dval, 1, queue );
        alpha = rho_new /ptq;
        
        
        magma_daxpy( dofs, alpha, p.dval, 1 , x->dval, 1, queue );                // x=x+alpha*p
        magma_daxpy( dofs, c_neg_one * alpha, q.dval, 1 , r.dval, 1, queue );     // r=r+alpha*q
        magma_daxpy( dofs, c_neg_one * MAGMA_D_CONJ(alpha), qt.dval, 1 , rt.dval, 1, queue );     // r=r+alpha*q

        res = magma_dnrm2( dofs, r.dval, 1, queue );

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_sync_wtime( queue );
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }

        if ( res/nomb <= solver_par->rtol || res <= solver_par->atol ){
            break;
        }
    }
    while ( solver_par->numiter+1 <= solver_par->maxiter );
    
    tempo2 = magma_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    double residual;
    CHECK(  magma_dresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->iter_res = res;
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter ) {
        info = MAGMA_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_SLOW_CONVERGENCE;
        if( solver_par->iter_res < solver_par->rtol*solver_par->init_res ||
            solver_par->iter_res < solver_par->atol ) {
            info = MAGMA_SUCCESS;
        }
    }
    else {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_DIVERGENCE;
    }
    
cleanup:
    magma_dmfree(&r, queue );
    magma_dmfree(&rt, queue );
    magma_dmfree(&p, queue );
    magma_dmfree(&pt, queue );
    magma_dmfree(&q, queue );
    magma_dmfree(&qt, queue );
    magma_dmfree(&y, queue );
    magma_dmfree(&yt, queue );
    magma_dmfree(&z, queue );
    magma_dmfree(&zt, queue );
    magma_dmfree(&AT, queue );
    magma_dmfree(&Ah1, queue );
    magma_dmfree(&Ah2, queue );

    solver_par->info = info;
    return info;
}   /* magma_dbicg */
