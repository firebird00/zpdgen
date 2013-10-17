#include <stdio.h>
#include <complex.h>
#include <gsl/gsl_vector.h>
#include "gpdf.h"

extern void wofzh_(double*, double*, double*, double*, int*);
extern void wofzwh_(double*, double*, double*, double*, int*);
extern void zbesj_(double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern void zbesi_(double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern void zbesk_(double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern void inmzpd_(double*, double*, double*,double *, int*, int*, int *, double*, double*, int*);
extern void inmkur_(double*, double*, double*,double *, int*, int*, double*, double*, int*);
extern void epszpd_(double*, double*, double*, int *, double*, double*, int*);
extern void zgeev_( char* jobvl, char* jobvr, int* n, complex* a,
		    int* lda, complex* w, complex* vl, int* ldvl, complex* vr, int* ldvr,
		    complex* work, int* lwork, double* rwork, int* info );

complex pldisp_inmzpd(complex za, double zb, double bi, int n, int m, int nw){
  complex res;
  double zr=creal(za);
  double zi=cimag(za);
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  int *flag=malloc(sizeof(int)*1);
  inmzpd_(&zr,&zi, &zb, &bi, &n,&m, &nw, cyr,cyi,flag);
  if(flag[0]!=0){
    printf("error in inmzpd!\n");
    printf("za=(%f,%f)\n",creal(za),cimag(za));
    exit(-1);
  }
  res=cyr[0]+I*cyi[0];
  free(cyr);
  free(cyi);
  return res;
}

complex pldisp_inmkur(complex za, double zb, double bi, int n, int m){
  complex res;
  double zr=creal(za);
  double zi=cimag(za);
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  int *flag=malloc(sizeof(int)*1);
  inmkur_(&zr,&zi, &zb, &bi, &n,&m, cyr,cyi,flag);
  if(flag[0]!=0){
    printf("error in inmzpd!\n");
    printf("za=(%f,%f)\n",creal(za),cimag(za));
    exit(-1);
  }
  res=cyr[0]+I*cyi[0];
  free(cyr);
  free(cyi);
  return res;
}

complex pldisp_epszpd(complex om, pldisp_eps_pars *pars,int nw){
  complex res;
  double zr=creal(om);
  double zi=cimag(om);
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  int *flag=malloc(sizeof(int)*1);
  epszpd_(&zr,&zi, (double *)pars, &nw, cyr,cyi,flag);
  if(flag[0]!=0){
    printf("error in epszpd!\n");
    printf("om=(%f,%f)\n",creal(om),cimag(om));
    exit(-1);
  }
  res=cyr[0]+I*cyi[0];
  free(cyr);
  free(cyi);
  return res;
}

complex pldisp_wofzh(complex z){
  complex res;
  double zr=creal(z);
  double zi=cimag(z);
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  int *flag=malloc(sizeof(int)*1);
  wofzh_(&zr,&zi,cyr,cyi,flag);
  if(flag[0]!=0){
    printf("error in wofz!\n");
    printf("z=(%f,%f)\n",creal(z),cimag(z));
    exit(-1);
  }
  res=cyr[0]+I*cyi[0];
  free(cyr);
  free(cyi);
  return res;
}

complex pldisp_wofzwh(complex z){
  complex res;
  double zr=creal(z);
  double zi=cimag(z);
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  int *flag=malloc(sizeof(int)*1);
  wofzwh_(&zr,&zi,cyr,cyi,flag);
  if(flag[0]!=0){
    printf("error in wofzw!\n");
    printf("z=(%f,%f)\n",creal(z),cimag(z));
    exit(-1);
  }
  res=cyr[0]+I*cyi[0];
  free(cyr);
  free(cyi);
  return res;
}

/*this is an interface to Amos' complex Bessel function implementation in fortran*/
complex pldisp_besselj(double nu, complex z){
  int kode=1;
  int n=1;
  double zr=creal(z);
  double zi=cimag(z);
  int nz;
  int *ierr=malloc(sizeof(int)*1);
  //  double cyr[1],cyi[1];
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  *cyr=0.0;
  *cyi=0.0;
  complex res;
  zbesj_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,ierr);
  if(ierr[0]!=0){
    printf("error in zbesj!\n");
    printf("z=(%f,%f)\n",creal(z),cimag(z));
    exit(-1);
  }
  res=cyr[0]+I*cyi[0];
  free(cyr);
  free(cyi);
  return res;
}

complex pldisp_besseli(double nu, complex z){
  int kode=1;
  int n=1;
  double zr=creal(z);
  double zi=cimag(z);
  int nz;
  int *ierr=malloc(sizeof(int)*1);
  //  double cyr[1],cyi[1];
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  *cyr=0.0;
  *cyi=0.0;
  complex res;
  zbesi_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,ierr);
  if(ierr[0]!=0){
    printf("error in zbesi!\n");
    exit(-1);
  }
  res=cyr[0]+I*cyi[0];
  free(cyr);
  free(cyi);
  return res;
}

complex pldisp_besselk(double nu, complex z){
  int kode=1;
  int n=1;
  double zr=creal(z);
  double zi=cimag(z);
  int nz;
  int *ierr=malloc(sizeof(int)*1);
  //  double cyr[1],cyi[1];
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  *cyr=0.0;
  *cyi=0.0;
  complex res;
  zbesk_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,ierr);
  if(ierr[0]!=0){
    printf("error in zbesk!\n");
  }
  res=cyr[0]+I*cyi[0];
  free(cyr);
  free(cyi);
  return res;
}


/* modified from numerical recipes' "gauher.c" -ODG, here is the original copyright note:*/
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
/* note that the indices run from 0 to n-1, in accord with the convention in c. 
   - Doesn't seem to work for anything larger than 128*/
/*
#define EPS_GH 3.0e-14
#define PIM4_GH 0.7511255444649425
#define MAXIT_GH 256
void pldisp_gauher(double *x,double *w, int n)
{
  int i,its,j,m;
  double p1,p2,p3,pp,z,z1;
  m=(n+1)/2;
  for (i=0;i<=m-1;i++) {
    if (i == 0) {
      z=sqrt((double)(2*n-1))-1.85575*pow((double)(2*n-1),-0.16667);
    } else if (i == 1) {
      z -= 1.14*pow((double)(n-1),0.426)/z;
    } else if (i == 2) {
      z=1.86*z-0.86*x[0];
    } else if (i == 3) {
      z=1.91*z-0.91*x[1];
    } else {
      z=2.0*z-x[i-2];
    }
    for (its=1;its<=MAXIT_GH;its++) {
      p1=PIM4_GH;
      p2=0.0;
      for (j=0;j<=n-1;j++) {
	p3=p2;
	p2=p1;
	p1=z*sqrt(2.0/(j+1.0))*p2-sqrt(((double)(j))/(j+1.0))*p3;
      }
      pp=sqrt((double)2*(n))*p2;
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= EPS_GH) break;
    }
    if (its > MAXIT_GH) printf("too many iterations in gauher\n");
    x[i]=z;
    x[n-1-i] = -z;
    w[i]=2.0/(pp*pp);
    w[n-1-i]=w[i];
  }
}
#undef EPS_GH
#undef PIM4_GH
#undef MAXIT_GH
*/
