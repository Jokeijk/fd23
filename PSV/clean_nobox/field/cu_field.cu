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
  safecall(cudaMalloc((void**)&(Txx),sizeof(float)*nx*tnz));
  safecall(cudaMalloc((void**)&(Txz),sizeof(float)*nx*tnz));
  safecall(cudaMalloc((void**)&(Tzz),sizeof(float)*nx*tnz));
  safecall(cudaMalloc((void**)&(U  ),sizeof(float)*nx*tnz));
  safecall(cudaMalloc((void**)&(W  ),sizeof(float)*nx*tnz));
  cudaMemset(Txx,  0,sizeof(float)*nx*tnz); Txx  -=(nz1-radius)*nx;
  cudaMemset(Txz,  0,sizeof(float)*nx*tnz); Txz  -=(nz1-radius)*nx;
  cudaMemset(Tzz,  0,sizeof(float)*nx*tnz); Tzz  -=(nz1-radius)*nx;
  cudaMemset(U  ,  0,sizeof(float)*nx*tnz); U    -=(nz1-radius)*nx;
  cudaMemset(W  ,  0,sizeof(float)*nx*tnz); W    -=(nz1-radius)*nx;
}

/* load from g_fld to gpu */
void FIELD::cu_load_restart(int deviceid,PARAM &param,FIELD & g_fld)
{
  nx=param.nx;
  nz=param.nz;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int tnz=nz2-nz1+1;
  safecall(cudaMemcpy(Txx+nz1*nx,g_fld.Txx+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(Txz+nz1*nx,g_fld.Txz+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(Tzz+nz1*nx,g_fld.Tzz+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(U  +nz1*nx,g_fld.U  +nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(W  +nz1*nx,g_fld.W  +nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
}

/* load from g_fld to gpu */
void FIELD::cu_save_state(int deviceid,PARAM &param,FIELD & g_fld)
{
  nx=param.nx;
  nz=param.nz;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int tnz=nz2-nz1+1;
  safecall(cudaMemcpy(g_fld.Txx+nz1*nx,Txx+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
  safecall(cudaMemcpy(g_fld.Txz+nz1*nx,Txz+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
  safecall(cudaMemcpy(g_fld.Tzz+nz1*nx,Tzz+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
  safecall(cudaMemcpy(g_fld.U  +nz1*nx,U  +nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
  safecall(cudaMemcpy(g_fld.W  +nz1*nx,W  +nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
}
