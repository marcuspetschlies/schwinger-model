#include "global.h"

void init_global ( void) {

L = 4;
T = 4;

V  = L*T;
Vh = V / 2;

xx_eo_up = NULL;
xx_eo_dn = NULL;

xx_lexic2eo = NULL;
xx_eo2lexic = NULL;

xx_ieo = NULL;

beta  = 1.;

m0 = 1.;

tau   = 1.;
nmd = 1;

deltaPhi = 0.5;

kappa = 0.;
kappa_bc = NULL;

verbose = 0;
}
