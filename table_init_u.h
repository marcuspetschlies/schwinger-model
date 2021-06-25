#ifndef _TABLE_INIT_U_H
#define _TABLE_INIT_U_H

/****************************************************
 * table_init_u.h
 *
 * PURPOSE:
 * DONE:
 * TODO:
 ****************************************************/

inline unsigned int * init_1level_utable ( size_t const N0 ) {
  return( N0 == 0 ? NULL : ( unsigned int *) calloc ( N0 , sizeof( unsigned int ) ) );
}  // end of init_1level_utable

/************************************************************************************/
/************************************************************************************/

inline void fini_1level_utable ( unsigned int **s  ) {
  if ( *s != NULL ) free ( *s );
  // fprintf ( stdout, "# [fini_1level_utable] active\n");
  *s = NULL;
}  // end of fini_1level_utable

/************************************************************************************/
/************************************************************************************/

inline unsigned int ** init_2level_utable (size_t const N0, size_t const N1 ) {
  unsigned int * s__ = NULL;
  s__ = init_1level_utable ( N0*N1);
  if ( s__ == NULL ) return( NULL );

  unsigned int ** s_ = ( N0 == 0 ) ? NULL : ( unsigned int **) malloc( N0 * sizeof( unsigned int *) );
  if ( s_ == NULL ) return ( NULL );

  for ( size_t i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_2level_utable

/************************************************************************************/
/************************************************************************************/


inline void fini_2level_utable ( unsigned int *** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_2level_utable] active\n");
    fini_1level_utable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_2level_utable

/************************************************************************************/
/************************************************************************************/


inline unsigned int *** init_3level_utable (size_t const N0, size_t const N1, size_t const N2 ) {
  unsigned int ** s__ = NULL;
  s__ = init_2level_utable ( N0*N1, N2);
  if ( s__ == NULL ) return( NULL );

  unsigned int *** s_ = ( N0 == 0 ) ? NULL : ( unsigned int ***) malloc( N0 * sizeof( unsigned int **) );
  if ( s_ == NULL ) return ( NULL );

  for ( size_t i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_3level_utable

/************************************************************************************/
/************************************************************************************/


inline void fini_3level_utable ( unsigned int **** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_3level_utable] active\n");
    fini_2level_utable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_3level_utable

/************************************************************************************/
/************************************************************************************/


inline unsigned int **** init_4level_utable (size_t const N0, size_t const N1, size_t const N2, size_t const N3 ) {
  unsigned int *** s__ = NULL;
  s__ = init_3level_utable ( N0*N1, N2, N3);
  if ( s__ == NULL ) return( NULL );

  unsigned int **** s_ = ( N0 == 0 ) ? NULL : ( unsigned int ****) malloc( N0 * sizeof( unsigned int ***) );
  if ( s_ == NULL ) return ( NULL );

  for ( size_t i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_4level_utable

/************************************************************************************/
/************************************************************************************/


inline void fini_4level_utable ( unsigned int ***** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_4level_utable] active\n");
    fini_3level_utable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_4level_utable

/************************************************************************************/
/************************************************************************************/


inline unsigned int ***** init_5level_utable (size_t const N0, size_t const N1, size_t const N2, size_t const N3, size_t const N4 ) {
  unsigned int **** s__ = NULL;
  s__ = init_4level_utable ( N0*N1, N2, N3, N4);
  if ( s__ == NULL ) return( NULL );

  unsigned int ***** s_ = ( N0 == 0 ) ? NULL : ( unsigned int *****) malloc( N0 * sizeof( unsigned int ****) );
  if ( s_ == NULL ) return ( NULL );

  for ( size_t i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_5level_utable

/************************************************************************************/
/************************************************************************************/


inline void fini_5level_utable ( unsigned int ****** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_5level_utable] active\n");
    fini_4level_utable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_5level_utable

/************************************************************************************/
/************************************************************************************/


inline unsigned int ****** init_6level_utable (size_t const N0, size_t const N1, size_t const N2, size_t const N3, size_t const N4, size_t const N5 ) {
  unsigned int ***** s__ = NULL;
  s__ = init_5level_utable ( N0*N1, N2, N3, N4, N5);
  if ( s__ == NULL ) return( NULL );

  unsigned int ****** s_ = ( N0 == 0 ) ? NULL : ( unsigned int ******) malloc( N0 * sizeof( unsigned int *****) );
  if ( s_ == NULL ) return ( NULL );

  for ( size_t i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_6level_utable

/************************************************************************************/
/************************************************************************************/


inline void fini_6level_utable ( unsigned int ******* s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_6level_utable] active\n");
    fini_5level_utable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_6level_utable

/************************************************************************************/
/************************************************************************************/


inline unsigned int ******* init_7level_utable (size_t const N0, size_t const N1, size_t const N2, size_t const N3, size_t const N4, size_t const N5, size_t const N6 ) {
  unsigned int ****** s__ = NULL;
  s__ = init_6level_utable ( N0*N1, N2, N3, N4, N5, N6);
  if ( s__ == NULL ) return( NULL );

  unsigned int ******* s_ = ( N0 == 0 ) ? NULL : ( unsigned int *******) malloc( N0 * sizeof( unsigned int ******) );
  if ( s_ == NULL ) return ( NULL );

  for ( size_t i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_7level_utable

/************************************************************************************/
/************************************************************************************/


inline void fini_7level_utable ( unsigned int ******** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_7level_utable] active\n");
    fini_6level_utable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_7level_utable

/************************************************************************************/
/************************************************************************************/


inline unsigned int ******** init_8level_utable (size_t const N0, size_t const N1, size_t const N2, size_t const N3, size_t const N4, size_t const N5, size_t const N6, size_t const N7 ) {
  unsigned int ******* s__ = NULL;
  s__ = init_7level_utable ( N0*N1, N2, N3, N4, N5, N6, N7);
  if ( s__ == NULL ) return( NULL );

  unsigned int ******** s_ = ( N0 == 0 ) ? NULL : ( unsigned int ********) malloc( N0 * sizeof( unsigned int *******) );
  if ( s_ == NULL ) return ( NULL );

  for ( size_t i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_8level_utable

/************************************************************************************/
/************************************************************************************/


inline void fini_8level_utable ( unsigned int ********* s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_8level_utable] active\n");
    fini_7level_utable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_8level_utable

/************************************************************************************/
/************************************************************************************/


#endif
