#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include "gpdf.h"
int main(int argc, char *argv[]){
  unsigned long int t,tc,dt;
  complex res;
  pldisp_hdf_file *outhdf,*outhdf2;
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
  dza=0.5+0.5i;
  omsi=-ps->ky;
  numx=creal(za_max-za_min)/creal(dza);
  numy=cimag(za_max-za_min)/cimag(dza);
  F=malloc(sizeof(complex)*numx*numy);
  za=malloc(sizeof(complex)*numx*numy);
  lt=0;
  for (k=0;k<4;k++){
    for(l=0;l<3;l++){
      n=nms[l*2];
      m=nms[l*2+1];
      printf("computing I%i%i for nw=%i\n",n,m,nws[k]);
      tc=clock();
      for(lx=0;lx<numx;lx++){
	for(ly=0;ly<numy;ly++){
	  za[lx*numy+ly]=za_min+creal(dza)*lx+I*cimag(dza)*ly;
	  F[lx*numy+ly]=pldisp_inmzpd(za[lx*numy+ly],0.0,0.3*0.3,n,m,nws[k]);
	}
      }
      dt=clock()-tc;
      dts[lt++]=dt/1.0e6;
      printf("\n%f secs cpu time\n",dts[lt-1]);
      sprintf(buf,"out_I%i%i_%i.h5",n,m,nws[k]);
      outhdf=pldisp_hdf5_create(buf);
      pldisp_hdf5_write_complex2d("za",outhdf,za,numx,numy);
      pldisp_hdf5_write_complex2d("eps",outhdf,F,numx,numy);
      H5Fclose(outhdf->id);
    }
  }

  for (k=0;k<4;k++){
    printf("computing eps for nw=%i\n",nws[k]);
    tc=clock();
    for(lx=0;lx<numx;lx++){
      for(ly=0;ly<numy;ly++){
	za[lx*numy+ly]=za_min+creal(dza)*lx+I*cimag(dza)*ly;
	om=-za[lx*numy+ly]*ps->omdi;
	I10=pldisp_inmzpd(za[lx*numy+ly],0.0,0.3*0.3,1,0,nws[k]);
	I30=pldisp_inmzpd(za[lx*numy+ly],0.0,0.3*0.3,3,0,nws[k]);
	I12=pldisp_inmzpd(za[lx*numy+ly],0.0,0.3*0.3,1,2,nws[k]);
	F[lx*numy+ly]=1.0+1.0/ps->tau+1.0/ps->omdi*((om-omsi*(1.0-1.5*ps->etai))*I10-omsi*ps->etai*(I30+I12));
      }
    }
    dt=clock()-tc;
    dts[lt++]=dt/1.0e6;
    printf("\n%f secs cpu time\n",dts[lt-1]);
    sprintf(buf,"out_eps_%i.h5",nws[k]);
    outhdf=pldisp_hdf5_create(buf);
    pldisp_hdf5_write_complex2d("za",outhdf,za,numx,numy);
    pldisp_hdf5_write_complex2d("eps",outhdf,F,numx,numy);
    H5Fclose(outhdf->id);
  }

  for (k=0;k<4;k++){
    printf("computing eps (combined) for nw=%i\n",nws[k]);
    tc=clock();
    for(lx=0;lx<numx;lx++){
      for(ly=0;ly<numy;ly++){
	za[lx*numy+ly]=za_min+creal(dza)*lx+I*cimag(dza)*ly;
	om=-za[lx*numy+ly]*ps->omdi;
	F[lx*numy+ly]=pldisp_epszpd(om,ps,nws[k]);
      }
    }
    dt=clock()-tc;
    dts[lt++]=dt/1.0e6;
    printf("\n%f secs cpu time\n",dts[lt-1]);
    sprintf(buf,"out_epsc_%i.h5",nws[k]);
    outhdf=pldisp_hdf5_create(buf);
    pldisp_hdf5_write_complex2d("za",outhdf,za,numx,numy);
    pldisp_hdf5_write_complex2d("eps",outhdf,F,numx,numy);
    H5Fclose(outhdf->id);
  }

  for(l=0;l<3;l++){
    n=nms[l*2];
    m=nms[l*2+1];
    printf("computing I%i%i (kuroda)\n",n,m);
    tc=clock();
    for(lx=0;lx<numx;lx++){
      for(ly=0;ly<numy;ly++){
	za[lx*numy+ly]=za_min+creal(dza)*lx+I*cimag(dza)*ly;
	F[lx*numy+ly]=pldisp_inmkur(za[lx*numy+ly],0.0,0.3*0.3,n,m);
      }
    }
    dt=clock()-tc;
    dts[lt++]=dt/1.0e6;
    printf("\n%f secs cpu time\n",dts[lt-1]);
    sprintf(buf,"out_I%i%i_kur.h5",n,m);
    outhdf=pldisp_hdf5_create(buf);
    pldisp_hdf5_write_complex2d("za",outhdf,za,numx,numy);
    pldisp_hdf5_write_complex2d("eps",outhdf,F,numx,numy);
    H5Fclose(outhdf->id);
  }
  printf("computing eps (kuroda)\n");
  tc=clock();
  for(lx=0;lx<numx;lx++){
    for(ly=0;ly<numy;ly++){
      za[lx*numy+ly]=za_min+creal(dza)*lx+I*cimag(dza)*ly;
      om=-za[lx*numy+ly]*ps->omdi;
      F[lx*numy+ly]=pldisp_epskur(om,ps);
    }
  }
  dt=clock()-tc;
  dts[lt++]=dt/1.0e6;
  printf("\n%f secs cpu time\n",dts[lt-1]);
  sprintf(buf,"out_epskur.h5");
  outhdf=pldisp_hdf5_create(buf);
  pldisp_hdf5_write_complex2d("za",outhdf,za,numx,numy);
  pldisp_hdf5_write_complex2d("eps",outhdf,F,numx,numy);
  H5Fclose(outhdf->id);
  outhdf2=pldisp_hdf5_create("time.h5");
  pldisp_hdf5_write_double1d("t",outhdf2,&dts[0],25);
  H5Fclose(outhdf2->id);
  //      F[lx*numy+ly]=I12;
  //0.5/ps->omdi*((om-omsi*(1.0-1.5*ps->etai))*I10-omsi*ps->etai*(I30+I12));
  //      F[lx*numy+ly]=pldisp_inmkur(za[lx*numy+ly],0.0,0.3*0.3,3,0);
  //            printf("za=%f %+fI\n",creal(za[lx*numy+ly]),cimag(za[lx*numy+ly]));
  //            printf("F=%f %+fI\n",creal(F[lx*numy+ly]),cimag(F[lx*numy+ly]));
  
  //    }
  //  }
}
