#ifndef _GLOBAL_H
#define _GLOBAL_H
#include <stdlib.h>
#include <stdio.h>

#if defined MAIN_PROGRAM
#  define EXTERN
#else
#  define EXTERN extern
#endif 

#define _GGI(_x,_k) (2*(_x)+(_k))
#define _GSI(_x) (4*(_x))

#ifndef CG_EPSREL
#define CG_EPSREL 1.e-16
#endif

#ifndef CG_MAXITR
#define CG_MAXITR 10000
#endif

/********************************************************************
 * global variables
 ********************************************************************/
EXTERN unsigned int L;
EXTERN unsigned int T;

EXTERN unsigned int V;
EXTERN unsigned int Vh;

EXTERN unsigned int ** xx_eo_up;
EXTERN unsigned int ** xx_eo_dn;

EXTERN unsigned int * xx_lexic2eo;
EXTERN unsigned int * xx_eo2lexic;

EXTERN unsigned int * xx_ieo;

EXTERN double beta;

EXTERN double m0;

EXTERN double tau;
EXTERN int nmd;

EXTERN double deltaPhi;

EXTERN double kappa;
EXTERN double ** kappa_bc;

EXTERN int verbose;


/***************************************************************************
 * calculate elapsed wall-time
 ***************************************************************************/
inline void show_time ( struct timeval * const ta, struct timeval * const tb, char * tag, char * timer, int const io ) {

  long int seconds =  tb->tv_sec  - ta->tv_sec;
  long int useconds = tb->tv_usec - ta->tv_usec;
  if ( useconds < 0 ) {
    useconds += 1000000;
    seconds--;
  }
  /* if ( io ) fprintf ( stdout, "# [%s] time for %s %ld sec %ld usec\n", tag, timer, seconds, useconds ); */
  if ( io ) fprintf ( stdout, "# [%s] time for %s %e sec\n", tag, timer, (double)seconds + (double)useconds/1000000.0 );

}  /* end of show_time */


#endif
