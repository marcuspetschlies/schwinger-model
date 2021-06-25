#ifndef _OPERATORS_H
#define _OPERATORS_H

void apply_D ( double * const y_out, double * const y_in, double * const g );

void apply_g5 ( double * const y , double * const x );

void apply_Ddag ( double * const y_out, double * const y_in, double * const g );

void apply_DDdag ( double * const y, double * const x, double * const g );

#endif
