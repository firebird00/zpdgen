///@file
#ifndef PLDISP_H
#define PLDISP_H
#include <complex.h>
#include <hdf5.h>

#define SQRT_TWO 1.41421356237310
#define SQRT_PI 1.77245385090552
#define EUL_GAM_TWO 0.2886078324507664327747
#define PLDISP_SHIFT_SINGULAR_LARGE 1e-4
#define PLDISP_SHIFT_SINGULAR_MEDIUM 1e-6
#define PLDISP_SHIFT_SINGULAR_SMALL 1e-8
#define PLDISP_LIMIT_SINGULAR 1e-8
#define PLDISP_LIMIT_SINGULAR_LARGE 1e-4
#define PLDISP_LIMIT_LARGEVAL 1e4
#define PLDISP_LARGEVAL 1e4

#define pldisp_signum(x) (((x)>=0.0) ? (1.0) : (-1.0))
#define pldisp_max(x,y) (((x)>(y)) ? (x) : (y))
#define pldisp_min(x,y) (((x)<(y)) ? (x) : (y))
#define pldisp_minus1_pow(m) ((m==((m/2)*2)) ? (1) : (-1))
#define pldisp_csqrt(z) (csqrt(-I*z)*(1.0+I)/SQRT_TWO)

#define PLDISP_MAXLISTELNUM 100

/// hdf file access structure
typedef struct pldisp_hdf_file_{
  hid_t id,gid;
}pldisp_hdf_file;

typedef struct pldisp_eps_pars_{
  double omdi,etai,tau;
  double ky,kpar;
}pldisp_eps_pars;

///the wrapper for the complex modified bessel function algorithm in fortran 77, due to Amos 83.
complex pldisp_besselj(double nu, complex z);
///the wrapper for the complex modified bessel function algorithm in fortran 77, due to Amos 83.
complex pldisp_besseli(double nu, complex z);
///the wrapper for the complex modified bessel function algorithm in fortran 77, due to Amos 83.
complex pldisp_besselk(double nu, complex z);
///the original wofz
complex pldisp_wofzh(complex z);
///the wofzw
complex pldisp_wofzwh(complex z);
///the inmzpd
complex pldisp_inmzpd(complex za, double zb, double bi, int n, int m, int nw);
///the inmkur
complex pldisp_inmkur(complex za, double zb, double bi, int n, int m);
///the epszpd
complex pldisp_epszpd(complex om,pldisp_eps_pars *ps, int nw);
/// utility function used to compute elapsed time.
unsigned long int pldisp_time();

pldisp_hdf_file * pldisp_hdf5_create(char *filename);
void pldisp_hdf5_write_double1d(char *varname, pldisp_hdf_file *fout, double *dat, int sz);
void pldisp_hdf5_write_complex1d(char *varname, pldisp_hdf_file *fout, complex *dat, int sz);
void pldisp_hdf5_write_double2d(char *varname, pldisp_hdf_file *fout, double *dat, int szx, int szy);
void pldisp_hdf5_write_complex2d(char *varname, pldisp_hdf_file *fout, complex *dat, int szx, int szy);

#endif
