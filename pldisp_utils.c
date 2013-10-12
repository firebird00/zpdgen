#include <sys/time.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <stddef.h>
#include <hdf5.h>
#include "gpdf.h"

unsigned long int pldisp_time(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec*1000+tv.tv_usec/1000;
}

pldisp_hdf_file * pldisp_hdf5_create(char *filename){
  pldisp_hdf_file *fout = malloc(sizeof(pldisp_hdf_file));
  fout->id=H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  fout->gid=H5Gcreate(fout->id,"fields",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  return fout;
}

void pldisp_hdf5_write_double1d(char *varname, pldisp_hdf_file *fout, double *dat, int sz){
  hid_t space_id,dtype_id,data_id,par_id;
  // space_id=H5Screate(H5S_SCALAR);
  hsize_t arsz[1];
  arsz[0]=sz;
  space_id=H5Screate_simple(1,arsz,arsz);
  par_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(par_id, 1, arsz);
  dtype_id=H5Tcopy(H5T_NATIVE_DOUBLE);
  //  H5Tset_size(dtype_id,arsz[0]);
  data_id=H5Dcreate(fout->gid,varname,dtype_id,space_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if(data_id<0){
    printf("could not create data object!");
  }
  H5Dextend(data_id,arsz);
  H5Dwrite(data_id,dtype_id,H5S_ALL, H5S_ALL, H5P_DEFAULT,dat);
  H5Dclose(data_id);
}

void pldisp_hdf5_write_double2d(char *varname, pldisp_hdf_file *fout, double *dat, int szx, int szy){
  hid_t space_id,dtype_id,data_id,par_id;
  // space_id=H5Screate(H5S_SCALAR);
  hsize_t arsz[2];
  arsz[0]=szx;
  arsz[1]=szy;
  space_id=H5Screate_simple(2,arsz,arsz);
  par_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(par_id, 2, arsz);
  dtype_id=H5Tcopy(H5T_NATIVE_DOUBLE);
  //  H5Tset_size(dtype_id,arsz[0]);
  data_id=H5Dcreate(fout->gid,varname,dtype_id,space_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if(data_id<0){
    printf("could not create data object!");
  }
  H5Dextend(data_id,arsz);
  H5Dwrite(data_id,dtype_id,H5S_ALL, H5S_ALL, H5P_DEFAULT,dat);
  H5Dclose(data_id);
}


void pldisp_hdf5_write_complex1d(char *varname, pldisp_hdf_file *fout, complex *dat, int sz){
  hid_t space_id,dtype_id,data_id,par_id;
  // space_id=H5Screate(H5S_SCALAR);
  hsize_t arsz[1];
  arsz[0]=sz;
  space_id=H5Screate_simple(1,arsz,arsz);
  par_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(par_id, 1, arsz);
  dtype_id = H5Tcreate (H5T_COMPOUND, sizeof (complex));
  H5Tinsert (dtype_id, "real", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert (dtype_id, "imag", sizeof(complex)/2, H5T_NATIVE_DOUBLE);
  data_id=H5Dcreate(fout->gid,varname,dtype_id,space_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if(data_id<0){
    printf("could not create data object!");
  }
  H5Dextend(data_id,arsz);
  H5Dwrite(data_id,dtype_id,H5S_ALL, H5S_ALL, H5P_DEFAULT,dat);
  H5Dclose(data_id);
}

void pldisp_hdf5_write_complex2d(char *varname, pldisp_hdf_file *fout, complex *dat, int szx,int szy){
  hid_t space_id,dtype_id,data_id,par_id;
  // space_id=H5Screate(H5S_SCALAR);
  hsize_t arsz[2];
  arsz[0]=szx;
  arsz[1]=szy;
  space_id=H5Screate_simple(2,arsz,arsz);
  par_id=H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(par_id, 2, arsz);
  dtype_id = H5Tcreate (H5T_COMPOUND, sizeof (complex));
  H5Tinsert (dtype_id, "real", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert (dtype_id, "imag", sizeof(complex)/2, H5T_NATIVE_DOUBLE);
  data_id=H5Dcreate(fout->gid,varname,dtype_id,space_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if(data_id<0){
    printf("could not create data object!");
  }
  H5Dextend(data_id,arsz);
  H5Dwrite(data_id,dtype_id,H5S_ALL, H5S_ALL, H5P_DEFAULT,dat);
  H5Dclose(data_id);
}
