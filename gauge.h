#ifndef _GAUGE_H
#define _GAUGE_H

double plaq ( double * const phi, int const dir );

double theta ( double * const phi, int const dir );

int read_config ( double * const c, char * const filename );

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

#endif
