#include"thread.h"
#include<cstdio>
#include<cstdlib>
#include"../param/param.h"
#include"../field/field.h"
#include"../gpu.h"
THREAD::THREAD(int ngpu)
{
  if(ngpu<=0) return;
  if(pthread_barrier_init(&barr, NULL, ngpu))
  {
	printf("Could not create a barrier\n");
	exit(-1);
  }
}
void THREAD::selectdevice(int deviceid)
{
  cudaSetDevice(deviceid);
}

void THREAD::wait(){
   pthread_barrier_wait(&barr);
}
inline int cpy_frombuf_boundary_ad( int deviceid,int nz1,int nz2,int nx,int nz,float*ad,float *gd)
{
  if(nz1>0){
    safecall(cudaMemcpy(ad+(nz1-radius)*nx,gd+(nz1-radius)*nx,sizeof(float)*nx*radius,cudaMemcpyHostToDevice));
  }
  if(nz2<nz-1){
    safecall(cudaMemcpy(ad+(nz2+1)*nx,gd+(nz2+1)*nx,sizeof(float)*nx*radius,cudaMemcpyHostToDevice));
  }
  return 0;
}
inline int cpy_backbuf_boundary_ad( int deviceid,int nz1,int nz2,int nx,int nz,float*ad,float *gd)
{
  if(nz1>0){
    safecall(cudaMemcpy(gd+(nz1)*nx,ad+(nz1)*nx,sizeof(float)*nx*radius,cudaMemcpyDeviceToHost));
  }
  if(nz2<nz-1){
    safecall(cudaMemcpy(gd+(nz2+1-radius)*nx,ad+(nz2+1-radius)*nx,sizeof(float)*nx*radius,cudaMemcpyDeviceToHost));
  }
  return 0;
}
int THREAD::exchange_stress(int deviceid,PARAM & param,FIELD & d_fld, FIELD &g_fld)
{
  if(param.ngpu==1) return 0;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int nx=g_fld.nx;
  int nz=g_fld.nz;
  cpy_backbuf_boundary_ad(deviceid,nz1,nz2,nx,nz,d_fld.Tx,g_fld.Tx);
  cpy_backbuf_boundary_ad(deviceid,nz1,nz2,nx,nz,d_fld.Tz,g_fld.Tz);
  this->wait();
  cpy_frombuf_boundary_ad(deviceid,nz1,nz2,nx,nz,d_fld.Tx,g_fld.Tx);
  cpy_frombuf_boundary_ad(deviceid,nz1,nz2,nx,nz,d_fld.Tz,g_fld.Tz);
  this->wait();
  return 0;
}
int THREAD::exchange_velocity(int deviceid,PARAM & param,FIELD & d_fld, FIELD &g_fld)
{
  if(param.ngpu==1) return 0;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int nx=g_fld.nx;
  int nz=g_fld.nz;
  cpy_backbuf_boundary_ad(deviceid,nz1,nz2,nx,nz,d_fld.V,g_fld.V);
  this->wait();
  cpy_frombuf_boundary_ad(deviceid,nz1,nz2,nx,nz,d_fld.V,g_fld.V);
  this->wait();
  return 0;
}
