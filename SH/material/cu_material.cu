#include"material.h"
#include"../param/param.h"
#include"../gpu.h"
#include<cstdio>

void MATERIAL::init_gpu_full(int deviceid,PARAM &param,MATERIAL &mat){
  nx=mat.nx;
  nz=mat.nz;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int tnz=nz2-nz1+1;

  safecall(cudaMemcpyToSymbol("d_coef",g_coef,sizeof(float)*20));

  usetable=mat.usetable;
  if(usetable){
	num_mat=mat.num_mat;
	usetable=mat.usetable;
	safecall(cudaMalloc((void**)&(tbl_BV ),sizeof(float)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_MUZ),sizeof(float)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_MUX),sizeof(float)*num_mat));
	safecall(cudaMalloc((void**)&(index),sizeof(float)*nx*tnz)); index  -=nz1*nx;

	safecall(cudaMemcpy( tbl_BV , mat.tbl_BV , sizeof(float)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_MUZ , mat.tbl_MUZ , sizeof(float)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_MUX, mat.tbl_MUX, sizeof(float)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( index  + nz1*nx, mat.index  + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));

  }else{
	safecall(cudaMalloc((void**)&(BV ),sizeof(float)*nx*tnz)); BV  -=nz1*nx;
	safecall(cudaMalloc((void**)&(MUX),sizeof(float)*nx*tnz));MUX  -=nz1*nx;
	safecall(cudaMalloc((void**)&(MUZ),sizeof(float)*nx*tnz)); MUZ -=nz1*nx;

	safecall(cudaMemcpy( BV  + nz1*nx, mat.BV  + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( MUX + nz1*nx, mat.MUX + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( MUZ + nz1*nx, mat.MUZ + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));

  }
}

