#ifndef _NR_UTILS_H
#define _NR_UTILS_H

#ifndef PI
#define PI acos(-1.0)
#endif
#define EPS 1.0e-7
#define MAX_ATOM 5000000
#define MAX_BP 100000

#define NR_END 1
#define FREE_ARG char*

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        	   (dmaxarg1) : (dmaxarg2))

static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        	   (dminarg1) : (dminarg2))

#include "nrutil_fncs.h"

#endif  /* _NR_UTILS_H */
