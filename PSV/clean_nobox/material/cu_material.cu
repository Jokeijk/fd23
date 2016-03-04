#include"material.h"
#include"../param/param.h"
#include"../gpu.h"
#include<cstdio>

extern __constant__ float d_coef[5][4];

void MATERIAL::init_gpu_full(int deviceid,PARAM &param,MATERIAL &mat){
  nx=mat.nx;
  nz=mat.nz;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int tnz=nz2-nz1+1;

  safecall(cudaMemcpyToSymbol(d_coef,g_coef,sizeof(float)*20));

  usetable=mat.usetable;
  if(usetable){
	num_mat=mat.num_mat;
	usetable=mat.usetable;
	safecall(cudaMalloc((void**)&(tbl_BU ),sizeof(float)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_BW ),sizeof(float)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_MU ),sizeof(float)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_MUA),sizeof(float)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_LAM),sizeof(float)*num_mat));
	safecall(cudaMalloc((void**)&(index),sizeof(float)*nx*tnz)); index  -=nz1*nx;

	safecall(cudaMemcpy( tbl_BU , mat.tbl_BU , sizeof(float)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_BW , mat.tbl_BW , sizeof(float)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_MU , mat.tbl_MU , sizeof(float)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_MUA, mat.tbl_MUA, sizeof(float)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_LAM, mat.tbl_LAM, sizeof(float)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( index  + nz1*nx, mat.index  + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));

  }else{
	safecall(cudaMalloc((void**)&(BU ),sizeof(float)*nx*tnz)); BU  -=nz1*nx;
	safecall(cudaMalloc((void**)&(BW ),sizeof(float)*nx*tnz)); BW  -=nz1*nx;
	safecall(cudaMalloc((void**)&(MU ),sizeof(float)*nx*tnz)); MU  -=nz1*nx;
	safecall(cudaMalloc((void**)&(MUA),sizeof(float)*nx*tnz)); MUA -=nz1*nx;
	safecall(cudaMalloc((void**)&(LAM),sizeof(float)*nx*tnz)); LAM -=nz1*nx;

	safecall(cudaMemcpy( BU  + nz1*nx, mat.BU  + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( BW  + nz1*nx, mat.BW  + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( MU  + nz1*nx, mat.MU  + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( MUA + nz1*nx, mat.MUA + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( LAM + nz1*nx, mat.LAM + nz1*nx, sizeof(float)*nx*tnz,cudaMemcpyHostToDevice));

  }

}
