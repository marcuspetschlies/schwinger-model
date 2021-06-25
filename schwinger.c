#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <complex.h>
#ifdef HAVE_MPI
#  include <mpi.h>
#endif
#ifdef HAVE_OPENMP
#  include <omp.h>
#endif
#include <getopt.h>

#ifndef TIMER
#define TIMER 1
#endif

#define MAIN_PROGRAM

#include "global.h"
#include "ranlxd.h"
#include "table_init_d.h"
#include "table_init_u.h"
#include "set_default.h"
#include "operators.h"
#include "solver.h"
#include "gauge.h"
#include "meas.h"


/********************************************************************/
/********************************************************************/

/*************************************************************************
 * pseudo random number generator for
 * Gaussian distribution
 * uses Box-Muller method
 * cf. Press, Teukolsky, Vetterling, Flannery: Numerical Receipes in C++. 
 * Second Edition. Cambridge University Press, 2002
 *************************************************************************/

int rangauss (double * const y1, unsigned int const NRAND) {

  double const TWO_MPI = 2. * M_PI;
  unsigned int const nrandh = NRAND/2;

  if(NRAND%2 != 0) {
    fprintf(stderr, "[rangauss] Error, NRAND must be an even number\n");
    return(1);
  }

  /* fill the complete field y1 */
  ranlxd(y1,NRAND);

#ifdef HAVE_OPEMP
#pragma omp parallel for 
#endif
  for ( unsigned int k = 0; k < nrandh; k++ ) {
    unsigned int const k2   = 2*k;
    unsigned int const k2p1 = k2+1;

    double const x1 = sqrt( -2. * log(y1[k2]) );
    y1[k2]   = x1 * cos( TWO_MPI * y1[k2p1] );
    y1[k2p1] = x1 * sin( TWO_MPI * y1[k2p1] );
  }  /* end of loop on nrandh */
  return(0);
}

/********************************************************************/
/********************************************************************/

/********************************************************************
 *
 ********************************************************************/
void spinor_rangauss ( double * const phi ) {
  double const norm = sqrt ( 0.5 );
  rangauss ( phi, _GSI(V) );
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < _GSI(V); ix++ ) {
    phi[ix] *= norm;
  }
}  /* end of spinor_rangauss */

/********************************************************************/
/********************************************************************/

void show_force ( double * const c, FILE * fs ) {
  int isopen = 0;
  if ( fs == NULL ) {
    fs = fopen ( "force.out", "w" );
    isopen = 1;
  }
  for ( unsigned int iy = 0; iy < V; iy++ ) {
    unsigned int const ix = xx_lexic2eo[iy];
    int const y0 = ( iy  / L );
    int const y1 = ( iy %  L );
    for ( int mu = 0; mu < 2; mu++) {
      fprintf ( fs, "%3d %3d     %d     %25.16e\n", y0, y1, mu, c[_GGI(ix, mu) ] );

    }
  }
  if( isopen ) fclose ( fs );

}  /* end of show_config */

/********************************************************************/
/********************************************************************/

void show_config ( double * const c, FILE * fs ) {
  int isopen = 0;
  if ( fs == NULL ) {
    fs = fopen ( "conf.out", "w" );
    isopen = 1;
  }
  for ( unsigned int iy = 0; iy < V; iy++ ) {
    unsigned int const ix = xx_lexic2eo[iy];
    int const y0 = ( iy  / L );
    int const y1 = ( iy %  L );
    for ( int mu = 0; mu < 2; mu++) {
      fprintf ( fs, "%3d %3d     %d     %25.16e %25.16e\n",
          y0, y1, mu, cos( c[_GGI(ix, mu) ] ), sin( c[_GGI(ix,mu) ] ) );

    }
  }
  if( isopen ) fclose ( fs );

}  /* end of show_config */

/********************************************************************/
/********************************************************************/

void show_fermion ( double * const c, FILE * fs ) {
  int isopen = 0;
  if ( fs == NULL ) {
    fs = fopen ( "fermion.out", "w" );
  }
  for ( unsigned int iy = 0; iy < V; iy++ ) {
    unsigned int const ix = xx_lexic2eo[iy];
    int const y0 = ( iy  / L );
    int const y1 = ( iy %  L );
    for ( int mu = 0; mu < 2; mu++) {
      fprintf ( fs, "%3d %3d     %d     %25.16e %25.16e\n", y0, y1, mu, c[_GSI(ix)+2*mu+0], c[_GSI(ix)+2*mu+1] );
    }
  }
  if( isopen ) fclose ( fs );

}  /* end of show_config */

/********************************************************************/
/********************************************************************/

/********************************************************************
 * complet Hamiltonian
 *
 * H = S_G + S_F
 *
 * S_G = beta x ( 1 - plaq ) x V
 * S_F = phi^+ ( M M^+ )^-1 phi
 *
 * plaq = volume-averaged plaquette
 * M = Dirac operator
 ********************************************************************/
double hamiltonian ( double * const phi, double * const p, double * const g ) {

  double h = beta * ( 1 - plaq ( g , 0) ) * V;

  h += hamiltonianf ( phi, p, g );

  return ( h );
}  /* end of hamiltonian */


/********************************************************************/
/********************************************************************/

void hmc_update ( double * const g ) {

  static double * p = NULL, * x = NULL, * phi = NULL;
  static unsigned int accept[2] = { 0, 0 };

  if ( g == NULL ) {
    /* reset */
    fini_1level_dtable ( &p );
    fini_1level_dtable ( &x );
    fini_1level_dtable ( &phi );
    fprintf ( stdout, "# [hmc_update] acceptance %6u %6u %16.7e\n", accept[0], accept[1], (double)(accept[0])/accept[1] );
    accept[0] = 0;
    accept[1] = 0;
    return;
  }

  if ( p == NULL ) {
    p   = init_1level_dtable ( 2*V );
    x   = init_1level_dtable ( 2*V );
    phi = init_1level_dtable ( _GSI(V) );
  }

  memcpy ( x, g, 2*V*sizeof(double) );

  rangauss ( p, 2*V );

  /* chi ~ exp ( - chî^+ chi ) */
  spinor_rangauss ( phi );

  /* phi = M chi 
   * in-place
   */
  apply_D ( phi, phi, g );

  /* Hamiltonian at start of trajectory */
  double const h0 = hamiltonian ( phi, p, x );

  /* MD leap frog evolution */
  leapfrog_update ( x, p, phi, tau, nmd );

  /* Hamiltonian at end of trajectory */
  double const h1 = hamiltonian ( phi, p, x );

  /* energy difference */
  double const dh = h1 - h0;

  double raccept;
  int acc=0;
  ranlxd ( &raccept, 1 );

  if ( raccept < exp( -dh ) ) {
    memcpy ( g, x, 2*V*sizeof(double) );
    acc = 1;
    accept[0]++;
  };
  accept[1]++;

  fprintf ( stdout, "hmc_update %16.7e %16.7e %d\n", h1, dh, acc );

  return ;
}

/********************************************************************/
/********************************************************************/


/********************************************************************
 * MAIN program
 ********************************************************************/
int main(int argc, char **argv) {

  int c;

  int Niter       = 0;
  int seed        = 0;
  int meas_every  = 1;
  double heat     = 0.;
  double bc_sign  = 0.;

  struct timeval ta, tb;

  init_global ();

  /********************************************************************
   * extract command line arguments
   ********************************************************************/
  while ((c = getopt(argc, argv, "hN:b:D:s:i:e:t:T:M:c:m:L:")) != -1) {
    switch (c) {
    case 'N':
      Niter = atoi ( optarg );
      fprintf ( stdout, "# [schwinger] Niter set to %d\n", Niter );
      break;
    case 's':
      seed = atoi ( optarg );
      fprintf ( stdout, "# [schwinger] seed set to %d\n", seed );
      break;
    case 'b':
      beta = atof ( optarg );
      fprintf ( stdout, "# [schwinger] beta set to %e\n", beta );
      break;
    case 'D':
      deltaPhi = atof ( optarg );
      fprintf ( stdout, "# [schwinger] deltaPhi set to %e\n", deltaPhi );
      break;
    case 'i':
      heat = atof ( optarg );
      fprintf ( stdout, "# [schwinger] heat set to %e\n", heat );
      break;
    case 'e':
      meas_every = atoi ( optarg );
      fprintf ( stdout, "# [schwinger] meas_every set to %d\n", meas_every );
      break;
    case 't':
      tau = atof ( optarg );
      fprintf ( stdout, "# [schwinger] tau set to %f\n", tau );
      break;
    case 'M':
      nmd = atoi ( optarg );
      fprintf ( stdout, "# [schwinger] nmd set to %d\n", nmd );
      break;
    case 'm':
      m0 = atof ( optarg );
      fprintf ( stdout, "# [schwinger] m0 set to %e\n", m0 );
      break;
    case 'c':
      bc_sign = atof ( optarg );
      fprintf ( stdout, "# [schwinger] bc sign set to %e\n", bc_sign );
      break;
    case 'T':
      T = atoi ( optarg );
      fprintf ( stdout, "# [schwinger] T set to %d\n", T );
      break;
    case 'L':
      L = atoi ( optarg );
      fprintf ( stdout, "# [schwinger] L set to %d\n", L );
      break;
    case 'h':
    case '?':
    default:
      fprintf ( stdout, "# [schwinger] Code for U1 sim\n" );
      exit (1);
      break;
    }
  }
  
  /********************************************************************
   * set volume
   ********************************************************************/
  V= T * L;
  Vh = V / 2;

  /********************************************************************
   * allocate geometry fields
   ********************************************************************/
  xx_eo_up = init_2level_utable ( V, 2 );
  xx_eo_dn = init_2level_utable ( V, 2 );

  xx_lexic2eo = init_1level_utable ( V );
  xx_eo2lexic = init_1level_utable ( V );
  
  xx_ieo = init_1level_utable ( V );

  unsigned int xeven = 0;
  unsigned int xodd  = 0;

  for ( unsigned int ix = 0; ix < V; ix++ ) {

    int const x0 = ix /  L;
    int const x1 = ix %  L;

    int const ieo = ( x0 + x1 ) % 2;

    if ( ieo == 0 ) {
      xx_lexic2eo[ix] = xeven;

      xx_eo2lexic[xeven] = ix;

      xx_ieo[ix] = 0;

      xeven++;
    } else {
      xx_lexic2eo[ix] = xodd + Vh;
      
      xx_eo2lexic[xodd+Vh] = ix;
      
      xx_ieo[ix] = 1;
      
      xodd++;
    }
  }

  /********************************************************************
   * nearest-neighbour fields
   ********************************************************************/
  for ( unsigned int ix = 0; ix < V; ix++ ) {

    int const x0 = ix / L;
    int const x1 = ix % L;

    /* pbc */
    int const y0 = ( x0 + 1 ) % T;
    int const y1 = ( x1 + 1 ) % L;

    int const z0 = ( x0 - 1 + T ) % T;
    int const z1 = ( x1 - 1 + L ) % L;

    xx_eo_up[ xx_lexic2eo[ix] ][0] = xx_lexic2eo[ y0 * L + x1 ];
    xx_eo_up[ xx_lexic2eo[ix] ][1] = xx_lexic2eo[ x0 * L + y1 ];

    xx_eo_dn[ xx_lexic2eo[ix] ][0] = xx_lexic2eo[ z0 * L + x1 ];
    xx_eo_dn[ xx_lexic2eo[ix] ][1] = xx_lexic2eo[ x0 * L + z1 ];

  }  /* end of loop on volume */

  /********************************************************************
   *
   ********************************************************************/
  kappa = 1. / ( 2. * m0 + 4. );
  kappa_bc = init_2level_dtable ( T , 2 );
  for ( unsigned int it = 0; it < T; it++ ) {
    kappa_bc[it][0] = kappa;
    kappa_bc[it][1] = kappa;
  }
  kappa_bc[  0][1] *= bc_sign;  /* timeslice   0 backward */
  kappa_bc[T-1][0] *= bc_sign;  /* timeslice T-1 forward  */

  /********************************************************************
   * initialize random number generator , Lüscher's RANLUX 
   * in double precision
   ********************************************************************/
  rlxd_init ( 2, seed);

  /********************************************************************
   * gauge configuration field
   *
   * = field of angles U = exp( I g )
   ********************************************************************/
  double * g = init_1level_dtable ( 2*V );

  /********************************************************************
   * heated start for g
   *
   * heat parameter interpolates between cold ( heat = 0 ) 
   * and fully random ( heat = 1. )
   ********************************************************************/
  ranlxd ( g, 2*V) ;
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++ ) {
    g[ix] = heat * 2. * M_PI * g[ix];
  }

  /* plaquette of initial gauge field */
  double plaq_g = plaq ( g , 0);
  fprintf ( stdout, "# [schwinger] plaq %4d %16.7e \n", 0, plaq_g );

  show_config ( g, NULL );

#if 0
  /********************************************************************
   * prepare random spinor field
   ********************************************************************/

  double ** spinor = init_2level_dtable ( 3, _GSI(V) );

  FILE * fs = NULL;

  ranlxd ( spinor[0], _GSI(V) ) ;
  // spinor[0][_GSI(xx_lexic2eo[0])] = 1.;

  fs = fopen ( "fermion.out", "w" );
  show_fermion ( spinor[0], fs );
  fclose ( fs );


  /********************************************************************
   * test apply_D
   ********************************************************************/
  apply_D ( spinor[1], spinor[0], g );

  fs = fopen ( "Dfermion.out", "w" );
  show_fermion ( spinor[1], fs );
  fclose ( fs );

  /********************************************************************
   * test apply_Ddag
   ********************************************************************/
  apply_Ddag ( spinor[2], spinor[0], g );

  fs = fopen ( "Ddagfermion.out", "w" );
  show_fermion ( spinor[2], fs );
  fclose ( fs );

  /********************************************************************
   * test apply_D
   ********************************************************************/
  apply_g5 ( spinor[1], spinor[0] );
  apply_D ( spinor[2], spinor[1], g );
  apply_g5 ( spinor[1], spinor[2] );

  fs = fopen ( "g5Dg5fermion.out", "w" );
  show_fermion ( spinor[1], fs );
  fclose ( fs );
#endif

#if 0
  memset ( spinor[1], 0, _GSI(V)*sizeof(double) );

  exitstatus = cg (spinor[1], spinor[0], g , CG_MAXITR, CG_EPSREL );

  fprintf ( stdout, "# [schwinger] exit status from cg was %d %s %d\n", exitstatus, __FILE__, __LINE__ );

  apply_Ddag( spinor[1], spinor[1], g );
  apply_D( spinor[1], spinor[1], g );

  fs = fopen ( "Difermion.out", "w" );
  show_fermion ( spinor[1], fs );
  fclose ( fs );

  fini_2level_dtable ( &spinor );
#endif

#if 0
  double * g_old = init_1level_dtable ( 2*V );
  memcpy ( g_old, g, 2*V*sizeof(double) );

  double * p     = init_1level_dtable ( 2*V );
  double * p_old = init_1level_dtable ( 2*V );

  rangauss ( p, 2*V);
  memcpy ( p_old, p, 2*V*sizeof(double) );

  double * chi = init_1level_dtable ( _GSI(V) );
  double * phi = init_1level_dtable ( _GSI(V) );

  /* chi ~ exp ( - chî^+ chi ) */
  spinor_rangauss ( chi );

  FILE * fs = fopen ( "chi.out", "w" );
  show_fermion ( chi, fs );
  fclose ( fs );

  /* phi = M chi */
  apply_D ( phi, chi, g );

  fs = fopen ( "phi.out", "w" );
  show_fermion ( phi, fs );
  fclose ( fs );

  double h_aux = 0;
  for ( unsigned int ix = 0; ix < _GSI(V); ix++ ) 
    h_aux += chi[ix] * chi[ix];

  for ( unsigned int ix = 0; ix < 2*V; ix++ ) 
    h_aux += p[ix] * p[ix] * 0.5;

  fprintf ( stdout, "# [schwinger] initial hamiltonian by hand %16.7e \n", h_aux );

  double h_old = hamiltonian ( phi, p_old, g_old );
  fprintf ( stdout, "# [schwinger] initial hamiltonian %16.7e \n", h_old );
#endif

#if 0
  double * force = init_1level_dtable ( 2*V );

  fermion_force ( force , phi, g_old );

  fs = fopen ( "force.out", "w" );
  show_force ( force, fs );
  fclose ( fs );

  double const epsilon = 1.e-8;

  g_old[_GGI(11,1)] += epsilon;

  h_old = - ( hamiltonian ( phi, p_old, g_old ) - h_old ) / epsilon;
  fprintf ( stdout, "# [schwinger] -d hamiltonian / d epsilon = %16.7e \n", h_old );

  fini_1level_dtable ( &force );

#endif

#if 0

  /********************************************************************
   * test integrator convergence
   ********************************************************************/

  tau = 1.;

  for ( nmd = 1; nmd <= (1<<16); nmd*=2 ) {

#if TIMER
    gettimeofday ( &ta, (struct timezone *)NULL );
#endif

    memcpy ( g, g_old, 2*V*sizeof(double) );

    memcpy ( p, p_old, 2*V*sizeof(double) );

    leapfrog_update ( g, p, phi, tau, nmd );

#if TIMER
    char timer_msg[200];
    sprintf ( timer_msg, "leapfrog_update-nmd=%d", nmd );
    gettimeofday ( &tb, (struct timezone *)NULL );
    show_time ( &ta, &tb, "schwinger", timer_msg, 1 );
#endif


    double h_new = hamiltonian ( phi, p, g );

    fprintf ( stdout, "ec %6.4f %4d %16.7e %16.7e %16.7e\n", tau, nmd, h_new, h_old, fabs ( ( h_new - h_old ) / h_old ) );

  }

  fini_1level_dtable ( &g_old );
  fini_1level_dtable ( &p );
  fini_1level_dtable ( &p_old );
  fini_1level_dtable ( &chi );
  fini_1level_dtable ( &phi );
#endif

  /********************************************************************
   *
   * Monte-Carlo iteration
   *
   ********************************************************************/
  FILE * meas_fs = fopen ( "meas.tab", "w" );

  int meas_count = 0;
  for ( int iter = 0; iter < Niter; iter++ ) {
    /* fprintf ( stdout, "# [u1] iter %d \n", iter ); */

#if TIMER
    gettimeofday ( &ta, (struct timezone *)NULL );
#endif
    /********************************************************************
     * HMC update
     ********************************************************************/
    hmc_update ( g );
#if TIMER
    gettimeofday ( &tb, (struct timezone *)NULL );
    show_time ( &ta, &tb, "u1", "update", 1 );
#endif

    /********************************************************************
     * measurement of observables
     ********************************************************************/
    if ( (iter +1) % meas_every == 0 ) {

      double const S[2] = { plaq ( g , 0), theta ( g, 0 ) };
      fprintf ( meas_fs, "%4d %4d %16.7e %16.7e\n", iter+1, meas_count, S[0], S[1] );

      meas_count++;

      pion_correlator ( g , 1, meas_fs );
    }
  }  /* end of loop on iterations */

#if 0
#endif
  leapfrog_update ( NULL, NULL, NULL, 0, 0 );
  hmc_update ( NULL );

  fclose ( meas_fs );


  fini_2level_dtable ( &kappa_bc );

  /********************************************************************
   * finalize
   ********************************************************************/

  fini_2level_utable ( &xx_eo_up );
  fini_2level_utable ( &xx_eo_dn );
  fini_1level_utable ( &xx_eo2lexic );
  fini_1level_utable ( &xx_lexic2eo );
  fini_1level_utable ( &xx_ieo );

  fini_1level_dtable ( &g );

  return ( 0 );
}
