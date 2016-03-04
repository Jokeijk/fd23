#include"field.h"
#include"../gpu.h"
#include"../param/param.h"
#include<cstdio>

/* init fld in CUDA */
void FIELD::init_gpu_full(int deviceid,PARAM & param){
  nx=param.nx;
  nz=param.nz;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int tnz=nz2-nz1+1+2*radius;
  safecall(cudaMalloc((void**)&(V  ),sizeof(float)*nx*tnz));
  safecall(cudaMalloc((void**)&(Tx ),sizeof(float)*nx*tnz));
  safecall(cudaMalloc((void**)&(Tz ),sizeof(float)*nx*tnz));
  cudaMemset(V  ,  0,sizeof(float)*nx*tnz); V    -=(nz1-radius)*nx;
  cudaMemset(Tx ,  0,sizeof(float)*nx*tnz); Tx   -=(nz1-radius)*nx;
  cudaMemset(Tz ,  0,sizeof(float)*nx*tnz); Tz   -=(nz1-radius)*nx;
}
