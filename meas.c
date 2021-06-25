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
#include "solver.h"

/********************************************************************/
/********************************************************************/

void spinor_ranz2 ( double * const r , unsigned int N) {
  double const norm = sqrt( 0.5 );

  ranlxd ( r, N );

#pragma omp parallel for
  for ( unsigned int k = 0; k < N; k++ ) {
      r[k] = (double)(2 * (int)(r[k]>=0.5) - 1) * norm;
  }
}  /* end of spinor_ranz2 */

/********************************************************************/
/********************************************************************/

void pion_correlator ( double * const g, int nsrc, FILE * fs ) {
  
  FILE * myfs = fs == NULL ? stdout : fs;

  double ** spinor = init_2level_dtable ( 2, _GSI(V) );
  double * corr = init_1level_dtable ( T );
  double * lran = init_1level_dtable ( _GSI(L) );

  for ( int isrc = 0; isrc < nsrc; isrc++ ) {

    memset ( spinor[0], 0, _GSI(V)*sizeof(double) );
    memset ( spinor[1], 0, _GSI(V)*sizeof(double) );

    double dtmp;
    ranlxd ( &dtmp, 1 );
    int const ts = (unsigned int)(T * dtmp ) % T;

    ranlxd ( lran, _GSI(L) );
    for ( unsigned int ix = 0; ix < L; ix++ ) {
      unsigned int const iy = _GSI( xx_lexic2eo[ts * L + ix] );
      memcpy ( spinor[0]+iy , lran+_GSI(ix) , _GSI(1)*sizeof(double) );
    }

    int exitstatus = cg ( spinor[1], spinor[0], g , CG_MAXITR, CG_EPSREL );

    for ( unsigned int it = 0; it < T; it++ ) {
      unsigned int itt = ( it + ts ) % T;
      double dtmp = 0.;
      for ( unsigned int ix = 0; ix < L; ix++ ) {
        unsigned int const iy = _GSI( xx_lexic2eo[itt * L + ix] );
        double * const _s = spinor[1]+iy;

        dtmp += _s[0] * _s[0] + _s[1] * _s[1] + _s[2] * _s[2] + _s[3] * _s[3];
      }
      corr[it] += dtmp;
    }
  }
  for ( unsigned int it = 0; it < T; it++ ) {
    fprintf ( myfs, "pion %3d %16.7e\n", it, corr[it] / (double)nsrc );
  }


  fini_2level_dtable ( &spinor );
  fini_1level_dtable ( &corr );
  fini_1level_dtable ( &lran );

}  /* end of pion_correlator */
