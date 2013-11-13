///@file
#ifndef GPDF_INM_H
#define GPDF_INM_H
#include <complex.h>

typedef struct gpdf_eps_pars_{
  double omdi,etai,tau;
  double ky,kpar;
}gpdf_eps_pars;

///c wrapper for inmzpd
complex gpdf_inm(complex za, double zb, double bi, int n, int m);
///c wrapper for epszpd
complex gpdf_eps(complex om, gpdf_eps_pars *pars);
#endif
