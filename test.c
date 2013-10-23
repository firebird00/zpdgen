#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include "gpdf.h"
int main(int argc, char *argv[]){
  unsigned long int t,tc,dt;
  complex res;
  pldisp_hdf_file *outhdf;
  int lx,ly;
  int n,m,lt,l,k;
  complex *F;
  complex *za,om,I10,I30,I12;
  complex za_min, za_max, dza;
  int numx,numy;
  pldisp_eps_pars *ps;
  double omsi;
  double dts[25];
  char buf[50];
  int nms[]={1,0,1,2,3,0};
  int nws[]={1,4,8,16};
  ps=malloc(sizeof(pldisp_eps_pars));
  ps->omdi=-0.5;
  ps->etai=3.4;
  ps->tau=1.0;
  ps->ky=0.3;
  ps->kpar=0.0;
  za_min=-6.0-6.05i;
  za_max=6.0+5.95i;
  dza=0.1+0.1i;
  omsi=-ps->ky;
  numx=creal(za_max-za_min)/creal(dza);
  numy=cimag(za_max-za_min)/cimag(dza);
  F=malloc(sizeof(complex)*numx*numy);
  za=malloc(sizeof(complex)*numx*numy);
  lt=0;
  for (k=0;k<4;k++){
    for(l=0;l<4;l++){
      n=nms[l*2];
      m=nms[l*2+1];
      printf("computing I%i%i for nw=%i\n",n,m,nws[k]);
      tc=clock();
      for(lx=0;lx<numx;lx++){
	for(ly=0;ly<numy;ly++){
	  za[lx*numy+ly]=za_min+creal(dza)*lx+I*cimag(dza)*ly;
	  om=za[lx*numy+ly];
	  F[lx*numy+ly]=pldisp_inmzpd(za[lx*numy+ly],0.0,0.3*0.3,n,m,nws[k]);
	}
      }
      dt=clock()-tc;
      dts[lt++]=dt/1.0e6;
      printf("\n%f secs cpu time\n",dts[lt-1]);
      sprintf(buf,"out_I%i%i.h5",n,m);
    }
  }
	  //      F[lx*numy+ly]=pldisp_epszpd(za[lx*numy+ly],ps);
	  //      I10=pldisp_inmzpd(om,0.0,0.3*0.3,1,0);
	  //      I30=pldisp_inmzpd(om,0.0,0.3*0.3,3,0);
	  //      I12=pldisp_inmzpd(om,0.0,0.3*0.3,1,2);
	  //      F[lx*numy+ly]=1.0+1.0/ps->tau+0.5/ps->omdi*((om-omsi*(1.0-1.5*ps->etai))*I10-omsi*ps->etai*(I30+I12));
	  //      F[lx*numy+ly]=I12;
	  //0.5/ps->omdi*((om-omsi*(1.0-1.5*ps->etai))*I10-omsi*ps->etai*(I30+I12));
	  //      F[lx*numy+ly]=pldisp_inmkur(za[lx*numy+ly],0.0,0.3*0.3,3,0);
	  //            printf("za=%f %+fI\n",creal(za[lx*numy+ly]),cimag(za[lx*numy+ly]));
	  //            printf("F=%f %+fI\n",creal(F[lx*numy+ly]),cimag(F[lx*numy+ly]));

      outhdf=pldisp_hdf5_create(buf);
      pldisp_hdf5_write_complex2d("za",outhdf,za,numx,numy);
      pldisp_hdf5_write_complex2d("eps",outhdf,F,numx,numy);
      H5Fclose(outhdf->id);
//    }
//  }
}
