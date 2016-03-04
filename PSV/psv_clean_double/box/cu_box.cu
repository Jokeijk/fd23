#include"box.h"
#include"../gpu.h"
#include"../field/field.h"
#include"../material/material.h"
#include"../param/param.h"
void BOX::init_gpu(int deviceid,BOX& box,PARAM &param)
{
  boxdevice=box.boxdevice;
  if(deviceid!=boxdevice) return; 
  bw=box.bw;
  start_x=box.start_x;
  start_z=box.start_z;
  n_add=box.n_add;
  n_sub=box.n_sub;


  safecall(cudaMalloc((void**)&(a_add_idx),sizeof(int)*n_add));
  safecall(cudaMalloc((void**)&(r_add_idx),sizeof(int)*n_add));
  safecall(cudaMalloc((void**)&(a_sub_idx),sizeof(int)*n_sub));
  safecall(cudaMalloc((void**)&(r_sub_idx),sizeof(int)*n_sub));

  safecall(cudaMemcpy(r_add_idx,box.r_add_idx,sizeof(int)*n_add,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(a_add_idx,box.a_add_idx,sizeof(int)*n_add,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(r_sub_idx,box.r_sub_idx,sizeof(int)*n_sub,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(a_sub_idx,box.a_sub_idx,sizeof(int)*n_sub,cudaMemcpyHostToDevice));

  p_sfield=new FIELD();
  p_sfield->init_gpu_box(deviceid,param);
  p_bmat  =new MATERIAL();
  p_bmat->init_gpu_box(deviceid,*(box.p_bmat));


}

void BOX::cu_copy_source(int deviceid,BOX &box)
{
  if(deviceid!=boxdevice) return;
  box.setboxsource();
  safecall(cudaMemcpy(p_sfield->Txx, box.p_sfield->Txx, sizeof(int)*bw*bw,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(p_sfield->Txz, box.p_sfield->Txz, sizeof(int)*bw*bw,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(p_sfield->Tzz, box.p_sfield->Tzz, sizeof(int)*bw*bw,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(p_sfield->U  , box.p_sfield->U  , sizeof(int)*bw*bw,cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(p_sfield->W  , box.p_sfield->W  , sizeof(int)*bw*bw,cudaMemcpyHostToDevice));
}



__global__ void ker_subsource(double *pa,double* pr,double* ps,int* a_sub_idx,int * r_sub_idx,int nsub)
{
  int ix = blockIdx.x*blockDim.x + threadIdx.x;
  if(ix<nsub){
	pr[r_sub_idx[ix]]=pa[a_sub_idx[ix]]-ps[r_sub_idx[ix]];
  }
}
void BOX::cu_subsource(int deviceid,double *pa,double* pr,double* ps)
{
  if(deviceid != boxdevice) return;
  int threadsPerBlock = n_sub;
  int blocksPerGrid = 1;
  ker_subsource<<<blocksPerGrid,threadsPerBlock>>>(pa,pr,ps,a_sub_idx,r_sub_idx,n_sub);
  CUT_CHECK_ERROR("call_cu_subsource");
}

__global__ void ker_addsource(double *pa,double* pr,double* ps,int* a_add_idx,int * r_add_idx,int nadd)
{
  int ix = blockIdx.x*blockDim.x + threadIdx.x;
  if(ix<nadd){
	pa[a_add_idx[ix]]=pr[r_add_idx[ix]]+ps[r_add_idx[ix]];
  }
}
void BOX::cu_addsource(int deviceid,double *pa,double* pr,double* ps)
{
  if(deviceid != boxdevice) return;
  int threadsPerBlock = n_add;
  int blocksPerGrid = 1;
  ker_addsource<<<blocksPerGrid,threadsPerBlock>>>(pa,pr,ps,a_add_idx,r_add_idx,n_add);
  CUT_CHECK_ERROR("call_cu_addsource");
}

