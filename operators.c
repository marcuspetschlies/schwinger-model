#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <complex.h>

#include "global.h"
#include "table_init_d.h"


/********************************************************************/
/********************************************************************/

void apply_D ( double * const y_out, double * const y_in, double * const g ) {

  double * const y = y_out;

  static double * x = NULL;

  if ( x == NULL ) {
    fprintf( stdout, "# [apply_D] allocate auxilliary x %s %d\n", __FILE__, __LINE__ );
    x = init_1level_dtable ( _GSI(V) );
  }

  memcpy ( x, y_in, _GSI(V) * sizeof ( double ) );

  int const gperm[2][4] = {
    { 2, 3, 0, 1 },
    { 3, 2, 1, 0 }
  };

  double const gsign[2][4] = {
    { 1.,  1.,  1., 1. }, 
    { 1., -1., -1., 1. }
  };

  double const norm = 1. / ( 2. * kappa );

#pragma omp parallel for
  for ( unsigned int iz = 0; iz < V; iz++ ) {
    int const it = xx_eo2lexic[iz] / L;

    double const hp[2][2] = { { kappa_bc[it][0], kappa_bc[it][1] }, { kappa, kappa } };

    /* fprintf( stdout, "# [apply_D] z = %6d t = %3d    norm = %16.7e hp = %16.7e     %16.7e    %16.7e    %16.7e\n", iz, it, norm, 
        hp[0][0], hp[0][1], hp[1][0], hp[1][1] ); */

    double * const _y = y + _GSI( iz );
    double * const _x = x + _GSI( iz );

    _y[0] = _x[0];
    _y[1] = _x[1];
    _y[2] = _x[2];
    _y[3] = _x[3];

    for ( int mu = 0; mu < 2; mu++ ) {

      unsigned int const izp = xx_eo_up[iz][mu];
      unsigned int const izm = xx_eo_dn[iz][mu];

      double const ufwd[2] = { cos ( g[_GGI(  iz, mu) ] ) * hp[mu][0],  sin ( g[_GGI( iz,mu) ] ) * hp[mu][0] }; 

      double const ubwd[2] = { cos ( g[_GGI( izm, mu) ] ) * hp[mu][1], -sin ( g[_GGI(izm,mu) ] ) * hp[mu][1] }; 

      double * const _xp = x + _GSI(izp);

      double const sfwd[2][2] = { 
          { _xp[0] - gsign[mu][0] * _xp[gperm[mu][0] ],
            _xp[1] - gsign[mu][1] * _xp[gperm[mu][1] ] },
          { _xp[2] - gsign[mu][2] * _xp[gperm[mu][2] ],
            _xp[3] - gsign[mu][3] * _xp[gperm[mu][3] ] } };

      double * const _xm = x + _GSI(izm);

      double const sbwd[2][2] = { 
          { _xm[0] + gsign[mu][0] * _xm[gperm[mu][0] ],
            _xm[1] + gsign[mu][1] * _xm[gperm[mu][1] ] },
          { _xm[2] + gsign[mu][2] * _xm[gperm[mu][2] ],
            _xm[3] + gsign[mu][3] * _xm[gperm[mu][3] ] } };
 
      /* per spin component and real, imaginary part */
      _y[0]  -=  ( sfwd[0][0] * ufwd[0] - sfwd[0][1] * ufwd[1] ) + ( sbwd[0][0] * ubwd[0] - sbwd[0][1] * ubwd[1] );
      _y[1]  -=  ( sfwd[0][0] * ufwd[1] + sfwd[0][1] * ufwd[0] ) + ( sbwd[0][0] * ubwd[1] + sbwd[0][1] * ubwd[0] );
      _y[2]  -=  ( sfwd[1][0] * ufwd[0] - sfwd[1][1] * ufwd[1] ) + ( sbwd[1][0] * ubwd[0] - sbwd[1][1] * ubwd[1] );
      _y[3]  -=  ( sfwd[1][0] * ufwd[1] + sfwd[1][1] * ufwd[0] ) + ( sbwd[1][0] * ubwd[1] + sbwd[1][1] * ubwd[0] );

    }

    _y[0] *= norm;
    _y[1] *= norm;
    _y[2] *= norm;
    _y[3] *= norm;

  }

}  /* end of apply_D */

/********************************************************************/
/********************************************************************/

/********************************************************************
 * g5 x spinor
 * g5 = tau3
 ********************************************************************/
void apply_g5 ( double * const y , double * const x ) {

#pragma omp parallel for
  for ( unsigned int iz = 0; iz < V; iz++ ) {
    double * const _y = y + _GSI(iz);
    double * const _x = x + _GSI(iz);

    _y[0] =   _x[0];
    _y[1] =   _x[1];
    _y[2] =  -_x[2];
    _y[3] =  -_x[3];
  }
}  /* end of apply_g5 */

/********************************************************************/
/********************************************************************/

/********************************************************************
 * M_nm = delta_nm - k sum_mu [ 
 *     (1 - sigma_mu) U_mu,n delta_n+mu,m
 *   + (1 + sigma_mu) U_mu,n-mu^+ delta_n-mu,m 
 *
 * (M)^+_nm = M_mn^+  = delta_nm - k sum_mu [ 
 *     (1 - sigma_mu) U_mu,m^+ delta_m+mu,n
 *   + (1 + sigma_mu) U_mu,m-mu delta_m-mu,n 
 * 
 * = delta_nm - k sum_mu [ 
 *     (1 - sigma_mu) U_mu,n-mu^+ delta_m+mu,n
 *   + (1 + sigma_mu) U_mu,n delta_m-mu,n 
 *
 ********************************************************************/
void apply_Ddag ( double * const y_out, double * const y_in, double * const g ) {

  double * const y = y_out;

  static double * x = NULL;

  if ( x == NULL ) {
    fprintf( stdout, "# [apply_Ddag] allocate auxilliary x %s %d\n", __FILE__, __LINE__ );
    x = init_1level_dtable ( _GSI(V) );
  }

  memcpy ( x, y_in, _GSI(V) * sizeof ( double ) );


  int const gperm[2][4] = {
    { 2, 3, 0, 1 },
    { 3, 2, 1, 0 }
  };

  double const gsign[2][4] = {
    { 1.,  1.,  1., 1. }, 
    { 1., -1., -1., 1. }
  };

  double const norm = 1. / ( 2. * kappa );

#pragma omp parallel for
  for ( unsigned int iz = 0; iz < V; iz++ ) {
    int const it = xx_eo2lexic[iz] / L;

    double const hp[2][2] = { { kappa_bc[it][0], kappa_bc[it][1] }, { kappa, kappa } };

    /* fprintf( stdout, "# [apply_D] t = %3d    norm = %e hp = %e     %e    %e    %e\n", it, norm, 
        hp[0][0], hp[0][1], hp[1][0], hp[1][1] ); */

    double * const _y = y + _GSI( iz );
    double * const _x = x + _GSI( iz );

    _y[0] = _x[0];
    _y[1] = _x[1];
    _y[2] = _x[2];
    _y[3] = _x[3];

    for ( int mu = 0; mu < 2; mu++ ) {

      unsigned int const izp = xx_eo_up[iz][mu];
      unsigned int const izm = xx_eo_dn[iz][mu];

      double const ufwd[2] = { cos ( g[_GGI(  iz, mu) ] ) * hp[mu][0],  sin ( g[_GGI( iz,mu) ] ) * hp[mu][0] }; 

      double const ubwd[2] = { cos ( g[_GGI( izm, mu) ] ) * hp[mu][1], -sin ( g[_GGI(izm,mu) ] ) * hp[mu][1] }; 

      double * const _xp = x + _GSI(izp);

      double const sfwd[2][2] = { 
          { _xp[0] + gsign[mu][0] * _xp[gperm[mu][0] ],
            _xp[1] + gsign[mu][1] * _xp[gperm[mu][1] ] },
          { _xp[2] + gsign[mu][2] * _xp[gperm[mu][2] ],
            _xp[3] + gsign[mu][3] * _xp[gperm[mu][3] ] } };

      double * const _xm = x + _GSI(izm);

      double const sbwd[2][2] = { 
          { _xm[0] - gsign[mu][0] * _xm[gperm[mu][0] ],
            _xm[1] - gsign[mu][1] * _xm[gperm[mu][1] ] },
          { _xm[2] - gsign[mu][2] * _xm[gperm[mu][2] ],
            _xm[3] - gsign[mu][3] * _xm[gperm[mu][3] ] } };
 
      /* per spin component and real, imaginary part */
      _y[0]  -=  ( sfwd[0][0] * ufwd[0] - sfwd[0][1] * ufwd[1] ) + ( sbwd[0][0] * ubwd[0] - sbwd[0][1] * ubwd[1] );
      _y[1]  -=  ( sfwd[0][0] * ufwd[1] + sfwd[0][1] * ufwd[0] ) + ( sbwd[0][0] * ubwd[1] + sbwd[0][1] * ubwd[0] );
      _y[2]  -=  ( sfwd[1][0] * ufwd[0] - sfwd[1][1] * ufwd[1] ) + ( sbwd[1][0] * ubwd[0] - sbwd[1][1] * ubwd[1] );
      _y[3]  -=  ( sfwd[1][0] * ufwd[1] + sfwd[1][1] * ufwd[0] ) + ( sbwd[1][0] * ubwd[1] + sbwd[1][1] * ubwd[0] );

    }

    _y[0] *= norm;
    _y[1] *= norm;
    _y[2] *= norm;
    _y[3] *= norm;

  }

}  /* end of apply_Ddag */

/********************************************************************/
/********************************************************************/

/********************************************************************
 * apply D after Ddag
 * Ddag on x into y
 * D in-place on y into y
 ********************************************************************/
void apply_DDdag ( double * const y, double * const x, double * const g ) {

  apply_Ddag ( y, x, g );
  apply_D ( y, y, g );

}  /* end of apply_DDdag */

/********************************************************************/
/********************************************************************/
