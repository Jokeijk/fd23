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
  safecall(cudaMalloc((void**)&(Txx),sizeof(double)*nx*tnz));
  safecall(cudaMalloc((void**)&(Txz),sizeof(double)*nx*tnz));
  safecall(cudaMalloc((void**)&(Tzz),sizeof(double)*nx*tnz));
  safecall(cudaMalloc((void**)&(U  ),sizeof(double)*nx*tnz));
  safecall(cudaMalloc((void**)&(W  ),sizeof(double)*nx*tnz));
  cudaMemset(Txx,  0,sizeof(double)*nx*tnz); Txx  -=(nz1-radius)*nx;
  cudaMemset(Txz,  0,sizeof(double)*nx*tnz); Txz  -=(nz1-radius)*nx;
  cudaMemset(Tzz,  0,sizeof(double)*nx*tnz); Tzz  -=(nz1-radius)*nx;
  cudaMemset(U  ,  0,sizeof(double)*nx*tnz); U    -=(nz1-radius)*nx;
  cudaMemset(W  ,  0,sizeof(double)*nx*tnz); W    -=(nz1-radius)*nx;
}

/* init from box fld in global to fld in CUDA */
void FIELD::init_gpu_box(int deviceid,PARAM &param){
  if(deviceid != param.boxdevice){
	return;
  }
  nx=param.bw;
  nz=param.bw;
  safecall(cudaMalloc((void**)&(Txx),sizeof(double)*nx*nz));
  safecall(cudaMalloc((void**)&(Txz),sizeof(double)*nx*nz));
  safecall(cudaMalloc((void**)&(Tzz),sizeof(double)*nx*nz));
  safecall(cudaMalloc((void**)&(U  ),sizeof(double)*nx*nz));
  safecall(cudaMalloc((void**)&(W  ),sizeof(double)*nx*nz));
  cudaMemset(Txx,  0,sizeof(double)*nx*nz); 
  cudaMemset(Txz,  0,sizeof(double)*nx*nz); 
  cudaMemset(Tzz,  0,sizeof(double)*nx*nz); 
  cudaMemset(U  ,  0,sizeof(double)*nx*nz); 
  cudaMemset(W  ,  0,sizeof(double)*nx*nz); 
}
