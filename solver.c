#include <stdlib.h>
#include <math.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "table_init_d.h"
#include "global.h"
#include "operators.h"
#include "gauge.h"

/*********************************************************************
 * y <- a x + b y
 *********************************************************************/
inline void axpby ( double const a , double * const x, double const  b, double * const y ) {
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < _GSI(V); ix++ ) {
    y[ix] = a * x[ix] + b * y[ix];
  }
}  /* end of axpby */

/*********************************************************************
 * x <- a x 
 *********************************************************************/
inline void ax ( double const a , double * const x ) {
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < _GSI(V); ix++ ) {
    x[ix] *= a;
  }
}  /* end of axpby */


/*********************************************************************
 * s <-x^+ y real part
 *********************************************************************/
inline double rdot ( double * const x, double * const y ) {

  double s = 0;

#ifdef HAVE_OPENMP
  omp_lock_t writelock;

  omp_init_lock(&writelock);

#pragma omp parallel
{
#endif

  double saccum = 0.;

#pragma omp for
  for ( unsigned int ix = 0; ix < _GSI(V); ix++ ) {
    saccum += x[ix] * y[ix];
  }

#ifdef HAVE_OPENMP
  omp_set_lock(&writelock);
#endif

  s += saccum;

#ifdef HAVE_OPENMP
  omp_unset_lock(&writelock);

}  /* end of parallel region */
#endif

  return (s);

}  /* end of rdot */

/*********************************************************************
 * cg iteration
 *********************************************************************/
int cg ( double * const x, double * const b, double * const g , int const maxiter, double const epsrel ) {

  double const DONE  =  1.;
  double const DMONE = -1.;
  double const DZERO =  0.;

  int status = 0 ;

  double ** spinor = init_2level_dtable ( 3, _GSI(V) );

  double * const r = spinor[0];
  double * const d = spinor[1];
  double * const z = spinor[2];

  double const norm_b = sqrt( rdot ( b, b)  );
  if ( verbose > 2 ) fprintf( stdout, "# [cg] norm of rhs = %e\n", norm_b );

  /* r = A x0 */
  apply_DDdag ( r, x, g );

  /* r <- b - r */
  axpby ( DONE, b, DMONE, r );

  /* d <- r */
  axpby ( DONE, r, DZERO, d );

  double res2 = rdot ( r, r );
  if ( verbose > 2 ) fprintf( stdout, "# [cg] initial res2 = %e %s %d\n", res2, __FILE__, __LINE__ );

  int iter = 0;

  for ( ; iter < maxiter && sqrt(res2) > epsrel * norm_b; iter++ ) {

    apply_DDdag ( z, d, g );

    double const alpha = res2 / rdot ( d, z );

    axpby ( alpha, d, DONE, x );

    axpby ( -alpha, z, DONE, r );

    double const rr = rdot ( r, r );

    double const beta = rr / res2;

    axpby ( DONE, r, beta, d );

    res2 = rr;

    if ( verbose > 1 ) fprintf ( stdout, "# [cg] iter %d res2 %e\n", iter, res2 );

  }

  if ( iter == maxiter && sqrt(res2) >= epsrel * norm_b ) {
    status = -1;
    fprintf ( stdout, "# [cg] CG did not converge\n" );
  } else {
    status = iter;
  }

  apply_DDdag ( r, x, g );
  axpby ( DMONE, b, DONE, r );

  double const res_true = sqrt( rdot ( r, r ) );

  fprintf ( stdout, "# [cg] true residue %16.7e %d\n", res_true, iter );

  fini_2level_dtable ( &spinor );

  return ( status );

}  /* end of cg */

/*********************************************************************
 * Hamiltonian
 *
 * IN:
 * phi = field ( e.g. M chi, where chi was Gaussian )
 * p   = conjugate momenta to gague field angles
 * g   = gauge field angles
 *********************************************************************/

double hamiltonianf ( double * const phi, double * const p, double * const g ) {
 
  double * psi = init_1level_dtable ( _GSI( V) );

  double h = 0;

  /* h = sum_x,mu p_x,mu^2 / 2 */
  for ( unsigned int ix = 0; ix < _GGI(V,0); ix++ ) {
    h += p[ix] * p[ix];
  }
  h *= 0.5;

  /* psi <- ( M Mdag)^-1 phi */
  cg ( psi, phi, g , CG_MAXITR, CG_EPSREL );

  /* h <- h + phi^+ psi = phi^+ ( M Mdag )^-1 phi */
  h += rdot ( phi, psi );

  fini_1level_dtable( &psi );

  return ( h );
  
}  /* end of hamiltonianf */

/********************************************************************
 * gauge force
 ********************************************************************/
void gauge_force  ( double * const f, double * const g ) {

  /********************************************************************
   *  S_G =  1 - Re ( U K ) = 1 - UK/2 - K^+U^+/2
   *
   *  fG = -dH/dw_x,mu = sum_l=t,x ( -Im ( U_x,mu K^l_x,mu )
   *     = -beta (   sin ( phi_x,mu + k^t_x,mu ) 
   *               + sin ( phi_x,mu + k^x_x,mu ) )
   ********************************************************************/
  double a[2] = { 0., 0. };
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < V; ix++ ) {
    for ( int mu = 0; mu < 2; mu++ ) {
      staples ( a, g, ix , mu );

      f[_GGI(ix,mu)] = -beta * ( sin ( g[_GGI(ix,mu) ] + a[0]) + sin ( g[_GGI(ix,mu)] + a[1] ) );
    }
  }

}  /* end of gauge force */

/********************************************************************/
/********************************************************************/


/*********************************************************************
 * fermion force
 *
 * phi is taken as a constant background field here
 *********************************************************************/
void fermion_force ( double * const f , double * const phi, double * const g ) {

  int const gperm[2][4] = {
    { 2, 3, 0, 1 },
    { 3, 2, 1, 0 }
  };

  double const gsign[2][4] = {
    { 1.,  1.,  1., 1. },
    { 1., -1., -1., 1. }
  };

  double ** spinor = init_2level_dtable ( 2, _GSI(V) );

  double * const psi = spinor[0];
  double * const xi  = spinor[1];

  /* psi = ( M Mdag)^-1 phi */
  cg ( psi, phi, g, CG_MAXITR, CG_EPSREL );

  /* xi = Mdag psi */
  apply_Ddag ( xi, psi, g );

#pragma omp parallel for
  for ( unsigned int iz = 0; iz < V; iz++ ) {

    int const it = xx_eo2lexic[iz] / L;

    /* kappa including boundary condition */
    double const hp[2][2] = { { kappa_bc[it][0], kappa_bc[it][1] }, { kappa, kappa } };

    for ( int mu = 0; mu < 2; mu++ ) {
      
      unsigned int izp = xx_eo_up[iz][mu];

      double const ufwd[2] = { cos ( g[_GGI(  iz, mu) ] ) * hp[mu][0],  sin ( g[_GGI( iz,mu) ] ) * hp[mu][0] };
      
      double const ubwd[2] = { ufwd[0], -ufwd[1] };  /* U^+ */
    
    
      double * const _xp = xi + _GSI(izp);

      double const sfwd[2][2] = {
          { _xp[0] - gsign[mu][0] * _xp[gperm[mu][0] ],
            _xp[1] - gsign[mu][1] * _xp[gperm[mu][1] ] },
          { _xp[2] - gsign[mu][2] * _xp[gperm[mu][2] ],
            _xp[3] - gsign[mu][3] * _xp[gperm[mu][3] ] } };

      double * const _xm = xi + _GSI(iz);

      double const sbwd[2][2] = {
          { _xm[0] + gsign[mu][0] * _xm[gperm[mu][0] ],
            _xm[1] + gsign[mu][1] * _xm[gperm[mu][1] ] },
          { _xm[2] + gsign[mu][2] * _xm[gperm[mu][2] ],
            _xm[3] + gsign[mu][3] * _xm[gperm[mu][3] ] } };

      double const rfwd[2][2] = {
        { ( sfwd[0][0] * ufwd[0] - sfwd[0][1] * ufwd[1] ) , ( sfwd[0][0] * ufwd[1] + sfwd[0][1] * ufwd[0] ) },
        { ( sfwd[1][0] * ufwd[0] - sfwd[1][1] * ufwd[1] ) , ( sfwd[1][0] * ufwd[1] + sfwd[1][1] * ufwd[0] ) } };
        
      double const rbwd[2][2] = {
        { ( sbwd[0][0] * ubwd[0] - sbwd[0][1] * ubwd[1] ) , ( sbwd[0][0] * ubwd[1] + sbwd[0][1] * ubwd[0] ) },
        { ( sbwd[1][0] * ubwd[0] - sbwd[1][1] * ubwd[1] ) , ( sbwd[1][0] * ubwd[1] + sbwd[1][1] * ubwd[0] ) } };

      double const pfwd[2][2] = {
        { psi[_GSI(iz)+0], psi[_GSI(iz)+1] },
        { psi[_GSI(iz)+2], psi[_GSI(iz)+3] } };
       
      double const pbwd[2][2] = {
        { psi[_GSI(izp)+0], psi[_GSI(izp)+1] },
        { psi[_GSI(izp)+2], psi[_GSI(izp)+3] } };
      
      /* force term from psi^+ [ ... ] xi 
       *
       * calculate only imaginar parts as
       * ( a + ib)^+ ( x +iy ) = -b x + a y
       */
      f[_GGI(iz,mu)] = 
        /* + psi^+_x [ (1 - sigma) U_x,mu xi_x+mu */
        +( -pfwd[0][1] * rfwd[0][0] + pfwd[0][0] * rfwd[0][1] 
           -pfwd[1][1] * rfwd[1][0] + pfwd[1][0] * rfwd[1][1] 
         ) \
        /* - psi^+_x+mu [ (1 + sigma) U_x,mu^+ xi_x */
        -( -pbwd[0][1] * rbwd[0][0] + pbwd[0][0] * rbwd[0][1] 
           -pbwd[1][1] * rbwd[1][0] + pbwd[1][0] * rbwd[1][1] 
         );

      /* normalize with 1 / kappa */
      f[_GGI(iz,mu)] /= kappa;
    
    }
  }

  fini_2level_dtable ( &spinor );

}  /* fermion_force */

/*********************************************************************/
/*********************************************************************/

/*********************************************************************
 *
 *********************************************************************/
void leapfrog_update ( double * const g , double * const p, double * const phi, double const tau, int const nmd ) {

  static double ** f = NULL;
  double const dt = tau / nmd;

  if ( nmd == 0 ) {
    fini_2level_dtable ( &f );
    return;
  }

  if ( f == NULL ) {
    f = init_2level_dtable ( 2, 2*V );
  }

  /* initial, half step P update */
  fermion_force ( f[0], phi, g );
  gauge_force   ( f[1],      g );

#pragma omp parallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++) {
    p[ix] += ( f[0][ix] + f[1][ix] ) * dt * 0.5;
  }

  /********************************************************************
   * loop on intermediate steps
   ********************************************************************/
  for ( int it = 0; it < nmd - 1; it++ ) {

    /* full step G update */
#pragma omp parallel for
    for ( unsigned int ix = 0; ix < 2*V; ix++) {
      g[ix] += p[ix] * dt;
    }

    /* full step P update */
    fermion_force ( f[0], phi, g );
    gauge_force   ( f[1],      g );

#pragma omp parallel for
    for ( unsigned int ix = 0; ix < 2*V; ix++) {
      p[ix] += ( f[0][ix] + f[1][ix] )* dt;
    }

  }  /* enf of loop on md steps */

  /* final, full step G update */
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++) {
    g[ix] += p[ix] * dt;
  }

  /* final, half step P final update */
  fermion_force ( f[0], phi, g ); 
  gauge_force   ( f[1],      g );

#pragma omp parallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++) {
    p[ix] += ( f[0][ix] + f[1][ix] ) * dt * 0.5;
  }

}  /* end of leapfrog_update */


/********************************************************************/
/********************************************************************/
