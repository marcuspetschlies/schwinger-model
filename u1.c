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

#define _GGI(_x,_k) (2*(_x)+(_k))


/********************************************************************
 * global variables
 ********************************************************************/
unsigned int const L = 8;
unsigned int const T = 8;

unsigned int const V  = L*T;
unsigned int const Vh = V / 2;

unsigned int ** xx_eo_up = NULL;
unsigned int ** xx_eo_dn = NULL;

unsigned int * xx_lexic2eo = NULL;
unsigned int * xx_eo2lexic = NULL;
  
unsigned int * xx_ieo = NULL;

double beta  = 1.;

double tau   = 1.;
int nmd = 1;

double deltaPhi = 0.5;

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

inline void staples ( double * const a, double * const phi, unsigned int ix , int const mu) {

  int const mup1 = ( mu + 1 ) % 2;

  unsigned int const ixp0 = xx_eo_up[ ix ][ mu ];
  /* unsigned int const ixm0 = xx_eo_dn[ ix ][ mu ]; */

  unsigned int const ixp1 = xx_eo_up[ ix ][ mup1 ];
  unsigned int const ixm1 = xx_eo_dn[ ix ][ mup1 ];

  unsigned int const ixpm1 = xx_eo_dn[ xx_eo_up[ ix ][ mu ] ][ mup1 ];

  /*
   *  staple pos mupK:
   *
   *  x + mu
   *      ____<_____  x + mu+mupK
   *   |            |
   *   |            |  
   *   ^            ^  S_K
   *   |            |
   *   |  ____>_____| 
   *  x               x + mupK
   *
   *  K = 1
   *
   */

  /* x --->  U( x, mu )  U( x + mu, mupK )       U( x+mupK, mu )^+     U( x, mupK )^+  ---> x */

  a[0] =                 phi[ _GGI(ixp0, mup1) ] - phi[ _GGI( ixp1, mu) ] - phi[ _GGI( ix, mup1) ];

  /*
   *  staple neg mupK
   *
   *  x -mupK + mu
   *            ____<_____   x + mu
   *           |              |
   *           |              |  
   *  S_(3+K)  \/             ^  
   *           |              |
   *           |____>_____    | 
   *        x -mupK          x 
   *
   *  K = 1
   *
   */

  /* x ---> U( x, mu )  U( x + mu - mupK , mupK )^+  U( x-mupK, mu )^+   U( x-mupK, mupK ) ---> x */

  a[1] =               -phi[ _GGI(ixpm1, mup1) ]     - phi[ _GGI( ixm1, mu) ] + phi[ _GGI( ixm1, mup1) ];

  return;

}  /* end of staples */

/********************************************************************/
/********************************************************************/

int mh_update ( double * const phi ) {

  static double * ran = NULL;
  
  static double * racc = NULL;

  int accepted = 0;

  if ( phi == NULL ) {
    if ( ran  != NULL ) free ( ran );
    if ( racc != NULL ) free ( racc );
    return ( 0 );
  }

  if ( ran == NULL ) {
    ran  = (double*)malloc( 2 * Vh * sizeof( double ) );
  }

  if ( racc == NULL ) {
    racc = (double*)malloc( 2 * Vh * sizeof( double ) );
  }

  /* even / odd sublattice */
  for ( int ieo = 0; ieo <= 1; ieo++ ) {
   
    unsigned int const offset = ieo * Vh;

    ranlxd ( ran,  2 * Vh );
    ranlxd ( racc, 2 * Vh );

    for ( int mu = 0; mu < 2; mu++ ) {
 
#pragma omp parallel for
      for ( unsigned int iy = 0; iy < Vh; iy++ ) {

        unsigned int const ix = iy + offset;

        double a[2] = { 0., 0. };

        staples ( a, phi, ix, mu );

        double const phi0 = phi[ _GGI( ix, mu) ];
        double const phi1 = deltaPhi * ( 2. * ran[ _GGI(iy, mu) ] - 1. ) + phi0;

        /* fprintf ( stdout, "# [u1] phi0 %e phi1 %e \n", phi0 , phi1 ); */

        double const dS = beta * (
          -cos( a[0] + phi1 ) + cos( a[0] + phi0 )
          -cos( a[1] + phi1 ) + cos( a[1] + phi0 ));

        /* accept or reject */
        if ( racc[ _GGI(iy, mu) ] < exp(-dS) ) {
          phi[ _GGI(ix, mu) ] = phi1;
          accepted++;
        }

      }  /* end of loop on directions */

    }  /* end of loop on sites */

  }  /* end of loop on even / odd sublattice */

  fprintf ( stdout, "# [u1] acceptance rate %e\n", (double)accepted / ( 2 * V ) );
  return (0);
}  /* end of update */


/********************************************************************/
/********************************************************************/

inline double plaquette_angle (double * const phi, unsigned int ix, int mu, int nu , int const dir) {

  unsigned int const ixpmu = xx_eo_up[ ix ][ mu ];
  unsigned int const ixpnu = xx_eo_up[ ix ][ nu ];

  unsigned int const ixmmu = xx_eo_dn[ ix ][ mu ];
  unsigned int const ixmnu = xx_eo_dn[ ix ][ nu ];

  unsigned int const ixpmumnu = xx_eo_dn[ ixpmu ][ nu ];

  unsigned int const ixmmupnu = xx_eo_dn[ ixpnu ][ mu ];

  unsigned int const ixmmumnu = xx_eo_dn[ ixmnu ][ mu ];


  /*   
   *   4 different geometries for plaquette calculation
   *
   *                   x + nu
   *                   0----<----0 x+mu+nu
   *                   |         |
   *                   |         |
   *                   |         |
   *                   0---->----0 
   *                   x         x+mu
   *
   *                   x         x + mu
   *                   0----<----0 
   *                   |         |
   *                   |         |
   *                   |         |
   *                   0---->----0 
   *                   x - nu    x + mu - nu
   *
   *    x - mu + nu   x + nu
   *        0----<----0 
   *        |         |
   *        |         ^
   *        |         |
   *        0---->----0 
   *     x - mu       x
   *
   *    x - mu        x
   *        0----<----0 
   *        |         |
   *        |         ^
   *        |         |
   *        0---->----0 
   *   x - mu - nu    x - nu
   *
   */

  double a = 0.;

  switch ( dir ) {
    case 0:
      a = ( phi[ _GGI(ix , mu ) ] + phi[ _GGI( ixpmu ,  nu ) ]    - phi[ _GGI( ixpnu , mu ) ] - phi[ _GGI( ix ,  nu ) ] );
      /* fprintf( stdout, "plaq %6d   %3d %3d     %d %d %25.16e    %25.16e\n",
          xx_eo2lexic[ix], 
          xx_eo2lexic[ix] / L, 
          xx_eo2lexic[ix] % L, 
          mu, nu, a, cos ( a ) );
          */
      return( a );
      break;
    case 1:
      return ( -phi[ _GGI(ix, mu )] + phi[ _GGI( ixpmumnu , nu ) ] + phi[ _GGI( ixmnu , mu ) ] - phi[ _GGI(ixmnu , nu ) ] );
      break;
    case 2:
      return(  phi[ _GGI(ix, nu ) ] - phi[ _GGI(ixmmupnu ,  mu ) ] - phi[ _GGI( ixmmu , nu) ] + phi[ _GGI( ixmmu , mu ) ] );
      break;
    case 3:
      return(  -phi[ _GGI(ixmmu , mu ) ] - phi[ _GGI(ixmmumnu, nu) ] + phi[ _GGI(ixmmumnu,  mu )] + phi[_GGI(ixmnu , nu ) ] );
      break;
  }

}  /* end of plaquette_angle */

/********************************************************************/
/********************************************************************/

double plaq ( double * const phi, int const dir ) {

  double S = 0.;

#ifdef HAVE_OPENMP
  omp_lock_t writelock;

  omp_init_lock(&writelock);

#pragma omp parallel
{
#endif

  double Saccum = 0;

#pragma omp for
  for ( unsigned int ix = 0; ix < V; ix++ ) {
     Saccum += cos( plaquette_angle ( phi, ix, 0, 1, dir ) );
  }

#ifdef HAVE_OPENMP
  omp_set_lock(&writelock);
#endif

  S += Saccum;

#ifdef HAVE_OPENMP
  omp_unset_lock(&writelock);

}  /* end of parallel region */

  omp_destroy_lock(&writelock);
#endif

  return ( S / ( (double)V ) );

}  /* end of plaq */

/********************************************************************/
/********************************************************************/

void wilson_loop ( double ** const wl, double * const phi, int const rmax, int const tmax ) {

  double * plaq = init_1level_dtable ( V );

#ifdef HAVE_OPENMP
  omp_lock_t writelock;

  omp_init_lock(&writelock);
#endif


  /* calculate plaquettes */
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < V; ix++ ) {
    /* e-o index */
    unsigned int const iy = xx_lexic2eo[ ix ];
    plaq[ix] = plaquette_angle ( phi, iy, 0, 1, 0 );
  }

  /* FILE * fs = fopen ( "plaq2", "w" );
  for ( unsigned int ix = 0; ix < V; ix++ ) {
    // unsigned int const iy = xx_lexic2eo[ ix ];
      fprintf( fs, "%6d  %3d %3d      %d %d   %25.16e %25.16e\n",
          ix, ix / L, ix % L ),
          0, 1, plaq[ix], cos ( plaq[ix] ) );
  }
  fclose ( fs );
  */

  /********************************************************************
   * Wilson loops
   *
   * as cumulative sums of plaquette angles
   *
   * = "tiling with plaquettes"
   *
   ********************************************************************/
#ifdef HAVE_OPENMP
#pragma omp parallel
{
#endif
  /* thread-local Wilson loop array for accumulation */
  double ** _wl = init_2level_dtable ( tmax, rmax );

#pragma omp for
  for ( unsigned int ix = 0; ix < V; ix++ ) 
  {
    
    /* ix shifted in t-direction, start at ix itself  */
    unsigned int ixp0 = ix;

    memset ( _wl[0], 0, tmax*rmax*sizeof(double) );

    for ( int it = 0; it < tmax; it++ ) 
    {
      /* initialize values for it to values from previous it */
      if ( it > 0 ) {
        memcpy ( _wl[it], _wl[it-1], rmax*sizeof(double) );
      }

      unsigned int ixp1 = ixp0;

      double _wr_accum = 0.;

      for ( int ir = 0; ir < rmax; ir++ ) {
      
        _wr_accum += plaq[ixp1];  /* t - x */

        /* add the plaquette at point shifted by 1 in direction k */
        _wl[it][ir] += _wr_accum;

        /* shift ixpk in k direction, go via eo coordinates */
        ixp1 = xx_eo2lexic[ xx_eo_up[ xx_lexic2eo[ixp1] ][1] ];

        /* fprintf( stdout, "ix %6u   %3d %3d it %d ir %d  %16.7e     %16.7e\n", 
            ix,
            ix /   L,
            ix %   L          ),
            it, ir,
            _wl[it][ir], 
           plaq[ixp1],);
        */
      }  /* end of loop on ir  */

      /* shift ixp0 in t-direction, go via eo coordinates */
      ixp0 = xx_eo2lexic[ xx_eo_up[ xx_lexic2eo[ixp0] ][0] ];

    }  /* end of loop on it */

#ifdef HAVE_OPENMP
    omp_set_lock(&writelock);
#endif

    for ( int k = 0; k < tmax*rmax ; k++ ) {
      wl[0][k] += cos( _wl[0][k] );
    }

#ifdef HAVE_OPENMP
    omp_unset_lock(&writelock);
#endif

  }  /* end of loop on ix */


  fini_2level_dtable ( &_wl );

#ifdef HAVE_OPENMP
}  /* end of parallel region */
#endif

#pragma omp parallel for
  for ( int k = 0; k < tmax*rmax ; k++ ) {
    wl[0][k] /= (double)V;
  }


  fini_1level_dtable ( &plaq );
}  /* end of wilson_loop */

/********************************************************************/
/********************************************************************/

/********************************************************************/
/********************************************************************/

int read_config ( double * const c, char * const filename ) {

  FILE * ifs = fopen ( filename, "r" );

  for ( unsigned int ix = 0; ix < V; ix++ ) {

    unsigned int const iy = xx_lexic2eo[ix];

    for ( int mu = 0; mu < 2; mu++ ) {
      double dtmp[2] = {0., 0.};

      // fscanf ( ifs, "%lf %lf\n", dtmp, dtmp+1 );
      fscanf ( ifs, "(%lf%lfj)\n", dtmp, dtmp+1 );

      c[4*iy+mu] = atan2 ( dtmp[1], dtmp[0] );

     /* fprintf ( stdout, "config read %6u %6u   %d    %25.16e  %25.16e %25.16e      %25.16e %25.16e\n", ix, iy, mu, 
         c[ _GGI(iy, mu) ], 
         cos( c[_GGI( iy, mu) ] ), 
         sin( c[ _GGI(iy, mu) ] ), 
         dtmp[0], dtmp[1] );
         */
    }
  }
  fclose ( ifs );
  
  return ( 0 );
}  /* end of read_config */


/********************************************************************/
/********************************************************************/

void show_config ( double * const c ) {
  for ( unsigned int iy = 0; iy < V; iy++ ) {
    unsigned int const ix = xx_lexic2eo[iy];
    int const y0 = ( iy  / L );
    int const y1 = ( iy %  L );
    for ( int mu = 0; mu < 2; mu++) {
      fprintf ( stdout, "config      %3d %3d     %d     %25.16e %25.16e\n",
          y0, y1, mu, cos( c[_GGI(ix, mu) ] ), sin( c[_GGI(ix,mu) ] ) );

    }
  }
}  /* end of show_config */

/********************************************************************/
/********************************************************************/

void gauge_force  ( double * const f, double * const phi ) {

/********************************************************************
 *
 *  S_G =  1 - Re ( U K ) = 1 - UK/2 - K^+U^+/2
 ********************************************************************/

  double a[2] = { 0., 0. };
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < V; ix++ ) {
    for ( int mu = 0; mu < 2; mu++ ) {
      staples ( a, phi, ix , mu );

      f[_GGI(ix,mu)] = beta * ( sin ( phi[_GGI(ix,mu) ] + a[0]) + sin ( phi[_GGI(ix,mu)] + a[1] ) );
    }
  }

}  /* end of gauge force */

/********************************************************************/
/********************************************************************/

void leapfrog_update ( double * const phi , double * const p, double const tau, int const nmd ) {

  static double * f = NULL;
  double const dt = tau / nmd;

  if ( nmd == 0 ) {
    fini_1level_dtable ( &f );
    return;
  }

  if ( f == NULL ) {
    f = init_1level_dtable ( 2*V );
  }

  /* P update, initial half-step */
  gauge_force ( f, phi );

#pragma omp parallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++) {
    p[ix] -= f[ix] * dt * 0.5;
  }

  /********************************************************************
   * loop on intermediate steps
   ********************************************************************/
  for ( int it = 0; it < nmd - 1; it++ ) {

    /* PHI update */
#pragma omp parallel for
    for ( unsigned int ix = 0; ix < 2*V; ix++) {
      phi[ix] += p[ix] * dt;
    }

    /* P update */
    gauge_force ( f, phi );

#pragma omp parallel for
    for ( unsigned int ix = 0; ix < 2*V; ix++) {
      p[ix] -= f[ix] * dt;
    }

  }  /* enf of loop on md steps */

  /* PHI final update */
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++) {
    phi[ix] += p[ix] * dt;
  }

  /* P final update, half-step */
  gauge_force ( f, phi );

#pragma omp parallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++) {
    p[ix] -= f[ix] * dt * 0.5;
  }

}  /* end of leapfrog_update */


/********************************************************************/
/********************************************************************/

double hamiltonian ( double * const phi, double * const p ) {
 
  double h = beta * ( 1 - plaq ( phi , 0) ) * V;
#pragma omp prallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++ ) {
    h += p[ix] * p[ix] * 0.5;
  }

  return ( h );
}
/********************************************************************/
/********************************************************************/

void hmc_update ( double * const phi ) { 
  static double * p = NULL, * x = NULL;

  if ( phi == NULL ) {
    fini_1level_dtable ( &p );
    fini_1level_dtable ( &x );
    return;
  }

  if ( p == NULL ) {
    p = init_1level_dtable ( 2*V );
    x = init_1level_dtable ( 2*V );
  }

  memcpy ( x, phi, 2*V*sizeof(double) );

  rangauss ( p, 2*V );

  double const h0 = hamiltonian ( x, p );

  leapfrog_update ( x, p, tau, nmd );

  double const h1 = hamiltonian ( x, p );

  double const dh = h1 - h0;

  double raccept;
  int acc=0; 
  ranlxd ( &raccept, 1 );

  if ( raccept < exp( -dh ) ) {
    memcpy ( phi, x, 2*V*sizeof(double) );
    acc = 1;
  };

  fprintf ( stdout, "# [hmc_update] %16.7e %16.7e %d\n", h1, dh, acc );

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
  int rmax        = L / 2 - 1;
  int tmax        = T / 2 - 1;

  int exitstatus = 0;
  struct timeval ta, tb;

  /********************************************************************
   * extract command line arguments
   ********************************************************************/
  while ((c = getopt(argc, argv, "hN:b:D:s:i:m:r:t:T:M:")) != -1) {
    switch (c) {
    case 'N':
      Niter = atoi ( optarg );
      fprintf ( stdout, "# [u1] Niter set to %d\n", Niter );
      break;
    case 's':
      seed = atoi ( optarg );
      fprintf ( stdout, "# [u1] seed set to %d\n", seed );
      break;
    case 'b':
      beta = atof ( optarg );
      fprintf ( stdout, "# [u1] beta set to %e\n", beta );
      break;
    case 'D':
      deltaPhi = atof ( optarg );
      fprintf ( stdout, "# [u1] deltaPhi set to %e\n", deltaPhi );
      break;
    case 'i':
      heat = atof ( optarg );
      fprintf ( stdout, "# [u1] heat set to %e\n", heat );
      break;
    case 'm':
      meas_every = atoi ( optarg );
      fprintf ( stdout, "# [u1] meas_every set to %d\n", meas_every );
      break;
    case 'r':
      rmax = atoi ( optarg );
      fprintf ( stdout, "# [u1] rmax set to %d\n", rmax );
      break;
    case 't':
      tmax = atoi ( optarg );
      fprintf ( stdout, "# [u1] tmax set to %d\n", tmax );
      break;
    case 'T':
      tau = atof ( optarg );
      fprintf ( stdout, "# [u1] tau set to %f\n", tau );
      break;
    case 'M':
      nmd = atoi ( optarg );
      fprintf ( stdout, "# [u1] nmd set to %d\n", nmd );
      break;
    case 'h':
    case '?':
    default:
      fprintf ( stdout, "# [u1] Code for U1 sim\n" );
      exit (1);
      break;
    }
  }

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
   * initialize random number generator , Lüscher's RANLUX 
   * in double precision
   ********************************************************************/
  rlxd_init ( 2, seed);


  /********************************************************************
   * gauge configuration field
   *
   * = field of angles U = exp( I phi )
   ********************************************************************/
  double * phi = init_1level_dtable ( 2*V );

  /********************************************************************
   * heated start for phi
   *
   * heat parameter interpolates between cold ( heat = 0 ) 
   * and fully random ( heat = 1. )
   ********************************************************************/
  ranlxd ( phi, 2*V) ;
#pragma omp parallel for
  for ( unsigned int ix = 0; ix < 2*V; ix++ ) {
    phi[ix] = heat * 2. * M_PI * phi[ix];
  }

  /* plaquette of initial gauge field */
  double plaq_phi = plaq ( phi , 0);
  fprintf ( stdout, "# [u1] plaq %4d %16.7e \n", 0, plaq_phi );


  /********************************************************************
   *
   * TEST leapfrog MD evolution error scaling
   *
   ********************************************************************/
#if 0
  double * phi_old = init_1level_dtable ( 2*V );
  memcpy ( phi_old, phi, 2*V*sizeof(double) );

  double * p     = init_1level_dtable ( 2*V );
  double * p_old = init_1level_dtable ( 2*V );

  ranlxd ( p, 2*V );

  memcpy ( p_old, p, 2*V*sizeof(double) );

  double h_old = hamiltonian (phi_old, p_old );

  tau = 1.;

  for ( nmd = 1; nmd <= (1<<12); nmd*=2 ) {

    memcpy ( phi, phi_old, 2*V*sizeof(double) );

    memcpy ( p, p_old, 2*V*sizeof(double) );

    leapfrog_update ( phi , p, tau, nmd );

    double h_new = hamiltonian ( phi, p );

    fprintf ( stdout, "%6.4f %4d %16.7e %16.7e %16.7e\n", tau, nmd, h_new, h_old, fabs ( ( h_new - h_old ) / h_old ) );

  }

  leapfrog_update ( NULL, NULL, 0, 0 );

  fini_1level_dtable ( &p );
  fini_1level_dtable ( &p_old );
  fini_1level_dtable ( &phi_old );

#endif

  double *** wl = init_3level_dtable ( Niter/meas_every, tmax, rmax );
  if ( wl == NULL ) {
    fprintf( stderr, "[u1] Error fron init_3level_dtable %s %d\n", __FILE__, __LINE__ );
    exit ( 1 );
  }

  /********************************************************************
   *
   * Monte-Carlo iteration
   *
   ********************************************************************/
  int meas_count = 0;
  for ( int iter = 0; iter < Niter; iter++ ) {
    /* fprintf ( stdout, "# [u1] iter %d \n", iter ); */
    
#if TIMER
    gettimeofday ( &ta, (struct timezone *)NULL );
#endif
    /********************************************************************
     * Metropolis-Hastings update
     ********************************************************************/
    // mh_update ( phi );
    hmc_update ( phi );
#if TIMER
    gettimeofday ( &tb, (struct timezone *)NULL );
    show_time ( &ta, &tb, "u1", "update", 1 );
#endif

    /********************************************************************
     * measurement of observables
     ********************************************************************/
    if ( (iter +1) % meas_every == 0 ) {

      double const S = plaq ( phi , 0);
      fprintf ( stdout, "plaq %4d %16.7e \n", iter+1 , S );

#if TIMER
      gettimeofday ( &ta, (struct timezone *)NULL );
#endif
      wilson_loop ( wl[meas_count], phi, rmax, tmax );
#if TIMER
      gettimeofday ( &tb, (struct timezone *)NULL );
      show_time ( &ta, &tb, "u1", "wilson_loop", 1 );
#endif
      meas_count++;
    }
  }  /* end of loop on iterations */

  /********************************************************************
   * write Wilson loop data 
   ********************************************************************/
  char filename[400];
  sprintf ( filename, "wilson-loop.L%d.T%d.beta%6.4f.dat", L, T, beta );
  FILE * wfs = fopen ( filename, "w" );
  if ( wfs == NULL ) {
    fprintf( stderr, "[u1] Error from fopen for file %s %s %d\n", filename, __FILE__, __LINE__ );
    exit ( 2 );
  }
  for ( int i = 0; i < meas_count; i++ ) {
    for ( int it = 0; it < tmax; it++ ) {
      for ( int ir = 0; ir < rmax; ir++ ) {
        fprintf ( wfs, "%4d %3d %3d %25.16e\n", i, it, ir, wl[i][it][ir] );
      }
    }
  }
  fclose ( wfs );
  fini_3level_dtable ( &wl );


  /********************************************************************
   * finalize
   ********************************************************************/

  mh_update ( NULL );

  fini_2level_utable ( &xx_eo_up );
  fini_2level_utable ( &xx_eo_dn );
  fini_1level_utable ( &xx_eo2lexic );
  fini_1level_utable ( &xx_lexic2eo );
  fini_1level_utable ( &xx_ieo );

  fini_1level_dtable ( &phi );

  return ( 0 );
}
