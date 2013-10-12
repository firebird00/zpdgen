#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include "gpdf.h"
int main(int argc, char *argv[]){
  unsigned long int t,tc;
  complex res;
  pldisp_hdf_file *outhdf;
  int lx,ly;
  int n,m;
  complex *F;
  complex *za,om,I10,I30,I12;
  complex za_min, za_max, dza;
  int numx,numy;
  pldisp_eps_pars *ps;
  double omsi;
  char buf[50];
  ps=malloc(sizeof(pldisp_eps_pars));
  ps->omdi=-0.5;
  ps->etai=3.4;
  ps->tau=1.0;
  ps->ky=0.3;
  ps->kpar=0.0;
  za_min=-6.0-6.0i;
  za_max=6.0+6.0i;
  dza=0.1+0.1i;
  omsi=-ps->ky;
  numx=creal(za_max-za_min)/creal(dza);
  numy=cimag(za_max-za_min)/cimag(dza);
  F=malloc(sizeof(complex)*numx*numy);
  za=malloc(sizeof(complex)*numx*numy);
  //  n=3;m=4;
  //  pldisp_inmzpd(4.0+4.0i,0.0,1.0,1,0);
  for(n=0;n<=4;n++){
    for(m=0;m<=4;m++){
      t=pldisp_time();
      tc=clock();
      for(lx=0;lx<numx;lx++){
	for(ly=0;ly<numy;ly++){
	  //  for(lx=0;lx<1;lx++){
	  //    for(ly=0;ly<1;ly++){
	  za[lx*numy+ly]=za_min+creal(dza)*lx+I*cimag(dza)*ly;
	  om=za[lx*numy+ly];
	  //      F[lx*numy+ly]=pldisp_epszpd(za[lx*numy+ly],ps);
	  F[lx*numy+ly]=pldisp_inmzpd(za[lx*numy+ly],0.0,0.3*0.3,n,m);
	  //      I10=pldisp_inmzpd(om,0.0,0.3*0.3,1,0);
	  //      I30=pldisp_inmzpd(om,0.0,0.3*0.3,3,0);
	  //      I12=pldisp_inmzpd(om,0.0,0.3*0.3,1,2);
	  //      F[lx*numy+ly]=1.0+1.0/ps->tau+0.5/ps->omdi*((om-omsi*(1.0-1.5*ps->etai))*I10-omsi*ps->etai*(I30+I12));
	  //      F[lx*numy+ly]=I12;
	  //0.5/ps->omdi*((om-omsi*(1.0-1.5*ps->etai))*I10-omsi*ps->etai*(I30+I12));
	  //      F[lx*numy+ly]=pldisp_inmkur(za[lx*numy+ly],0.0,0.3*0.3,3,0);
	  //            printf("za=%f %+fI\n",creal(za[lx*numy+ly]),cimag(za[lx*numy+ly]));
	  //            printf("F=%f %+fI\n",creal(F[lx*numy+ly]),cimag(F[lx*numy+ly]));
	}
      }
      printf("\n%lu msecs real time\n",pldisp_time()-t);
      printf("\n%lu msecs cpu time\n",clock()-tc);
      sprintf(buf,"out_I%i%i.h5",n,m);
      outhdf=pldisp_hdf5_create(buf);
      pldisp_hdf5_write_complex2d("za",outhdf,za,numx,numy);
      pldisp_hdf5_write_complex2d("eps",outhdf,F,numx,numy);
      H5Fclose(outhdf->id);
    }
  }
}
