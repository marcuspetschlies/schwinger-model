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
#include "operators.h"
#include "solver.h"

/*********************************************************************
 * fill spinor with Z2xZ2 iid
 *********************************************************************/
void spinor_ranz2 ( double * const r , unsigned int N) {
  double const norm = sqrt( 0.5 );

  ranlxd ( r, N );

#pragma omp parallel for
  for ( unsigned int k = 0; k < N; k++ ) {
      r[k] = (double)(2 * (int)(r[k]>=0.5) - 1) * norm;
  }
}  /* end of spinor_ranz2 */

/*********************************************************************
 * calculate pseudoscalar / pseudovector correlators
 *********************************************************************/
void pion_correlator ( double * const g, int nsrc, FILE * fs ) {
  
  FILE * myfs = fs == NULL ? stdout : fs;

  double ** spinor = init_2level_dtable ( 2, _GSI(V) );
  double ** corr = init_2level_dtable (2, T );
  double * lran = init_1level_dtable ( _GSI(L) );

  /* loop on stochastic samples */
  for ( int isrc = 0; isrc < nsrc; isrc++ ) {

    /* init all to zero fields */
    memset ( spinor[0], 0, _GSI(V)*sizeof(double) );
    memset ( spinor[1], 0, _GSI(V)*sizeof(double) );

    /*********************************************************************
     * select random source timeslice
     *********************************************************************/
    double dtmp;
    ranlxd ( &dtmp, 1 );
    int const ts = (unsigned int)(T * dtmp ) % T;

    /*********************************************************************
     * fill timeslice with a random vector
     *********************************************************************/
    ranlxd ( lran, _GSI(L) );
#pragma omp parallel for
    for ( unsigned int ix = 0; ix < L; ix++ ) {
      unsigned int const iy = _GSI( xx_lexic2eo[ts * L + ix] );
      memcpy ( spinor[0]+iy , lran+_GSI(ix) , _GSI(1)*sizeof(double) );
    }

    /*********************************************************************
     * run inverter on random source, produce
     * spinor1 = (M M⁺)^-1 spinor0
     *********************************************************************/
    int exitstatus = cg ( spinor[1], spinor[0], g , CG_MAXITR, CG_EPSREL );
    if ( exitstatus < 0 ) {
      fprintf ( stderr, "[pion_correlator] Error from cg, status was %d %s %d\n", exitstatus, __FILE__, __LINE__);
      exit(1);
    }
    /* check status ... */

    /*********************************************************************
     * apply M^+ to get
     * spinor0 = M^+ spinor1 = M^+ (M M^+)^-1 spinor0 = M^-1 spinor0
     *********************************************************************/
    apply_Ddag ( spinor[0], spinor[1], g );

    /*********************************************************************
     * contract to correlators
     *********************************************************************/
#pragma omp parallel for
    for ( unsigned int it = 0; it < T; it++ ) {
      unsigned int itt = ( it + ts ) % T;
      double dtmp = 0., dtmp2 = 0.;
      for ( unsigned int ix = 0; ix < L; ix++ ) {
        unsigned int const iy = _GSI( xx_lexic2eo[itt * L + ix] );
        double * const _s = spinor[0]+iy;

        /* pion correlator, phi^+ phi for all x */
        dtmp  += ( _s[0] * _s[0] + _s[1] * _s[1] ) + ( _s[2] * _s[2] + _s[3] * _s[3] );

        /* a0 pi correlator, phi^+ g0 phi; g0 = ( 0 & 1 \\ 1 & 0 ) */
        dtmp2 += ( _s[0] * _s[2] + _s[1] * _s[3] ) + ( _s[2] * _s[0] + _s[3] * _s[1] );
      }

      corr[0][it] += dtmp;
      corr[1][it] += dtmp2;
    }
  }
  for ( unsigned int it = 0; it < T; it++ ) {
    fprintf ( myfs, "%3d %16.7e %16.7e\n", it, corr[0][it] / (double)nsrc, corr[1][it] / (double)nsrc );
  }


  fini_2level_dtable ( &spinor );
  fini_2level_dtable ( &corr );
  fini_1level_dtable ( &lran );

}  /* end of pion_correlator */
