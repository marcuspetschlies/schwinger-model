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

#include "global.h"
#include "ranlxd.h"
#include "table_init_d.h"
#include "table_init_u.h"

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

  return ( 0 );
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

/********************************************************************/
/********************************************************************/

int read_config ( double * const c, char * const filename ) {

  FILE * ifs = fopen ( filename, "r" );

  for ( unsigned int ix = 0; ix < V; ix++ ) {

    unsigned int const iy = xx_lexic2eo[ix];

    for ( int mu = 0; mu < 2; mu++ ) {
      double dtmp[2] = {0., 0.};

      // fscanf ( ifs, "%lf %lf\n", dtmp, dtmp+1 );
      if ( fscanf ( ifs, "(%lf%lfj)\n", dtmp, dtmp+1 ) != 2 ) {
        fprintf ( stderr, "[read_config] Error from fscanf %s %d\n", __FILE__, __LINE__ );
        return ( 1 );
      }

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

double theta ( double * const phi, int const dir ) {

  double const TWO_MPI = 2. * M_PI;
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
    double const dtmp  = plaquette_angle ( phi, ix, 0, 1, dir );
    double const dtmp2 = dtmp - ( (int)(dtmp / TWO_MPI) ) * TWO_MPI;
    double const dtmp3 = dtmp2 >= 0 ? dtmp2 : dtmp2 + TWO_MPI;
    double const dtmp4 = dtmp3 <= M_PI ? dtmp3 : dtmp3 - TWO_MPI;
    /* fprintf ( stdout, "phi-mod %16.7e %16.7e %16.7e %16.7e\n", dtmp, dtmp2, dtmp3, dtmp4); */
    Saccum += dtmp4;
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

  return ( S / TWO_MPI );

}  /* end of plaq */
