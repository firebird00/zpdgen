#include <stdio.h>
#include <complex.h>
#include "gpdf.h"

extern void inmzpd_(double*, double*, double*,double *, int*, int*, double*, double*, int*);

complex gpdf_inm(complex za, double zb, double bi, int n, int m){
  complex res;
  double zr=creal(za);
  double zi=cimag(za);
  double *cyr=malloc(sizeof(double)*1);
  double *cyi=malloc(sizeof(double)*1);
  int *flag=malloc(sizeof(int)*1);
  inmzpd_(&zr,&zi, &zb, &bi, &n,&m, cyr,cyi,flag);
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
