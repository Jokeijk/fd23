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

  safecall(cudaMemcpyToSymbol("d_coef",g_coef,sizeof(double)*20));

  usetable=mat.usetable;
  if(usetable){
	num_mat=mat.num_mat;
	usetable=mat.usetable;
	safecall(cudaMalloc((void**)&(tbl_BU ),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_BW ),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_MU ),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_MUA),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_LAM),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(index),sizeof(double)*nx*tnz)); index  -=nz1*nx;

	safecall(cudaMemcpy( tbl_BU , mat.tbl_BU , sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_BW , mat.tbl_BW , sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_MU , mat.tbl_MU , sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_MUA, mat.tbl_MUA, sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_LAM, mat.tbl_LAM, sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( index  + nz1*nx, mat.index  + nz1*nx, sizeof(double)*nx*tnz,cudaMemcpyHostToDevice));

  }else{
	safecall(cudaMalloc((void**)&(BU ),sizeof(double)*nx*tnz)); BU  -=nz1*nx;
	safecall(cudaMalloc((void**)&(BW ),sizeof(double)*nx*tnz)); BW  -=nz1*nx;
	safecall(cudaMalloc((void**)&(MU ),sizeof(double)*nx*tnz)); MU  -=nz1*nx;
	safecall(cudaMalloc((void**)&(MUA),sizeof(double)*nx*tnz)); MUA -=nz1*nx;
	safecall(cudaMalloc((void**)&(LAM),sizeof(double)*nx*tnz)); LAM -=nz1*nx;

	safecall(cudaMemcpy( BU  + nz1*nx, mat.BU  + nz1*nx, sizeof(double)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( BW  + nz1*nx, mat.BW  + nz1*nx, sizeof(double)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( MU  + nz1*nx, mat.MU  + nz1*nx, sizeof(double)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( MUA + nz1*nx, mat.MUA + nz1*nx, sizeof(double)*nx*tnz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( LAM + nz1*nx, mat.LAM + nz1*nx, sizeof(double)*nx*tnz,cudaMemcpyHostToDevice));

  }

}

void MATERIAL::init_gpu_box(int deviceid,MATERIAL &mat){
  nx=mat.nx;
  nz=mat.nz;
  usetable=mat.usetable;
  if(usetable){
	num_mat=mat.num_mat;
	usetable=mat.usetable;
	safecall(cudaMalloc((void**)&(tbl_BU ),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_BW ),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_MU ),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_MUA),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(tbl_LAM),sizeof(double)*num_mat));
	safecall(cudaMalloc((void**)&(index),sizeof(double)*nx*nz));

	safecall(cudaMemcpy( tbl_BU , mat.tbl_BU , sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_BW , mat.tbl_BW , sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_MU , mat.tbl_MU , sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_MUA, mat.tbl_MUA, sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( tbl_LAM, mat.tbl_LAM, sizeof(double)*num_mat,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( index  , mat.index  , sizeof(double)*nx*nz,cudaMemcpyHostToDevice));
  }else{
	safecall(cudaMalloc((void**)&(BU ),sizeof(double)*nx*nz));
	safecall(cudaMalloc((void**)&(BW ),sizeof(double)*nx*nz));
	safecall(cudaMalloc((void**)&(MU ),sizeof(double)*nx*nz));
	safecall(cudaMalloc((void**)&(MUA),sizeof(double)*nx*nz));
	safecall(cudaMalloc((void**)&(LAM),sizeof(double)*nx*nz));

	safecall(cudaMemcpy( BU  , mat.BU  , sizeof(double)*nx*nz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( BW  , mat.BW  , sizeof(double)*nx*nz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( MU  , mat.MU  , sizeof(double)*nx*nz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( MUA , mat.MUA , sizeof(double)*nx*nz,cudaMemcpyHostToDevice));
	safecall(cudaMemcpy( LAM , mat.LAM , sizeof(double)*nx*nz,cudaMemcpyHostToDevice));
  }
}
