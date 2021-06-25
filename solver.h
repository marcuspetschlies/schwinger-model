#ifndef _SOLVER_H
#define _SOLVER_H

int cg ( double * const x, double * const b, double * const g , int const maxiter, double const epsrel );

double hamiltonianf ( double * const phi, double * const p, double * const g );

void fermion_force ( double * const f , double * const phi, double * const g );

void leapfrog_update ( double * const g , double * const p, double * const chi, double const tau, int const nmd );

#endif
