#include"fd.h"
#include"../gpu.h"

#include"../param/const.h"
#include"../param/param.h"
#include"../material/material.h"
#include"../field/field.h"
#include"../cpml/cpml.h"
#include"../slide/slide.h"

#include<cstdio>
#include<cstdlib>

//#define topabsorb false

#define s_data(i,j) s_data[j][i]

#define CHECK_IF_NEED_UPDATE  \
  if(usepunch){\
	if(\
		current_step <= first_arrive_step[(iz/BDIMY)*(nx/BDIMX)+(ix/BDIMX)] ||\
		current_step >= last_step        [(iz/BDIMY)*(nx/BDIMX)+(ix/BDIMX)] ){ \
	  return;\
	}\
  }


#define  TxDIFF \
  Tx_x=d_coef[ord][0]*( s_data(tx,ty)   - s_data(tx-1,ty) )\
    -d_coef[ord][1]*( s_data(tx+1,ty) - s_data(tx-2,ty) )\
    +d_coef[ord][2]*( s_data(tx+2,ty) - s_data(tx-3,ty) )\
    -d_coef[ord][3]*( s_data(tx+3,ty) - s_data(tx-4,ty) );

#define TzDIFF \
  Tz_z=d_coef[ord][0]*( s_data(tx,ty+1) - s_data(tx,ty) )\
    -d_coef[ord][1]*( s_data(tx,ty+2) - s_data(tx,ty-1) )\
    +d_coef[ord][2]*( s_data(tx,ty+3) - s_data(tx,ty-2) )\
    -d_coef[ord][3]*( s_data(tx,ty+4) - s_data(tx,ty-3) );

#define  DxV \
  dxV=d_coef[ord][0]*( s_data(tx+1,ty)   - s_data(tx,ty) )\
       -d_coef[ord][1]*( s_data(tx+2,ty) - s_data(tx-1,ty) )\
       +d_coef[ord][2]*( s_data(tx+3,ty) - s_data(tx-2,ty) )\
       -d_coef[ord][3]*( s_data(tx+4,ty) - s_data(tx-3,ty) );


#define  DzV \
  dzV=d_coef[ord][0]*( s_data(tx,ty)   - s_data(tx,ty-1) )\
       -d_coef[ord][1]*( s_data(tx,ty+1) - s_data(tx,ty-2) )\
       +d_coef[ord][2]*( s_data(tx,ty+2) - s_data(tx,ty-3) )\
       -d_coef[ord][3]*( s_data(tx,ty+3) - s_data(tx,ty-4) );

#define ZERO_HOLO \
s_data[threadIdx.y][threadIdx.x]=0.0;\
s_data[threadIdx.y][BDIMX+2*radius-1-threadIdx.x]=0.0;\
s_data[BDIMY+2*radius-1-threadIdx.y][threadIdx.x]=0.0;\
s_data[BDIMY+2*radius-1-threadIdx.y][BDIMX+2*radius-1-threadIdx.x]=0.0; 

#define EDGE_SHARE(d_F) \
if(threadIdx.y<radius && iz>radius){\
	s_data[threadIdx.y][tx] = d_F[in_idx-radius*nx];\
  }\
if(threadIdx.y<radius && iz+BDIMY<nz){\
  s_data[threadIdx.y+BDIMY+radius][tx] = d_F[in_idx+BDIMY*nx];\
}\
if(threadIdx.x<radius && ix>radius){\
  s_data[ty][threadIdx.x] = d_F[in_idx-radius];\
}\
if(threadIdx.x<radius && ix+BDIMX<nx){\
  s_data[ty][threadIdx.x+BDIMX+radius] = d_F[in_idx+BDIMX];\
}\
s_data[ty][tx] = d_F[in_idx];

__constant__ float d_coef[5][4];


__global__ void cudaupdate_stress(
	bool usetable,
	int * userestrict mat_index,
	bool usepunch,
	int current_step,
	int * userestrict first_arrive_step,
	int * userestrict last_step,
	float * userestrict Tx,
	float * userestrict Tz,
	float * userestrict V,
	float * userestrict MUX,
	float * userestrict MUZ,
	float * userestrict psi_V_x,
	float * userestrict psi_V_z,
	float * userestrict b_V_x,
    float * userestrict b_V_z,
	float * userestrict c_V_x,
    float * userestrict c_V_z,
	float * userestrict k_V_x,
    float * userestrict k_V_z,
	int nx,
	int nz,
	int npml,
	int izshift)
{
  __shared__ float s_data[BDIMY+2*radius][BDIMX+2*radius];
  int ix = blockIdx.x*blockDim.x + threadIdx.x;
  int iz = blockIdx.y*blockDim.y + threadIdx.y + izshift;

  CHECK_IF_NEED_UPDATE;

  int ord;
  int idx;
  int in_idx = iz*nx + ix; // index for global memory
  int tx=threadIdx.x+radius;  //index for shared memory
  int ty=threadIdx.y+radius;
  float dxV,dzV;
  float MUX_ix_iz;
  float MUZ_ix_iz;

  if(usetable){
	MUX_ix_iz= MUX[mat_index[in_idx]];
	MUZ_ix_iz= MUZ[mat_index[in_idx]];
  }else{
	MUX_ix_iz= MUX[in_idx];
	MUZ_ix_iz= MUZ[in_idx];
  }

  ZERO_HOLO;
  __syncthreads();

    
  ord=min(ix,iz);
//  ord=ix; // maybe top layer can also use 8th,this reduced dispersion
  ord=min(ord,nx-1-ix);
  ord=min(ord,nz-1-iz);
  ord=min(ord,ORD8);

  EDGE_SHARE(V);
  __syncthreads();
  DxV; 
  DzV; 

  if(ix==0 || ix>=nx-1 || iz>=nz-1) return;

  if(ix<npml){
    //left pml
    idx=ix+iz*2*npml;
    psi_V_x[idx] = b_V_x[ix] * psi_V_x[idx] + c_V_x[ix] * dxV;
    Tx[in_idx] += MUX_ix_iz*( dxV * k_V_x[ix] + psi_V_x[idx] );

  }else if(ix>=nx-npml){
    //right pml
    idx=npml+nx-1-ix+iz*2*npml;
    psi_V_x[idx] = b_V_x[ix] * psi_V_x[idx] + c_V_x[ix] * dxV;
    Tx[in_idx] += MUX_ix_iz*( dxV * k_V_x[ix] + psi_V_x[idx] );
  }else{
    Tx[in_idx] += MUX_ix_iz*( dxV );
  }

  if(iz<npml && topabsorb){
	//top pml
	idx=(iz*nx)+ix;
    psi_V_z[idx] = b_V_z[iz] * psi_V_z[idx] + c_V_z[iz] * dzV;
    Tz[in_idx] += MUZ_ix_iz*( dzV * k_V_z[iz] + psi_V_z[idx] ) ;
  }else if(iz>=nz-npml){
	// bottom
	idx=(npml+nz-1-iz)*nx+ix;
    psi_V_z[idx] = b_V_z[iz] * psi_V_z[idx] + c_V_z[iz] * dzV;
    Tz[in_idx] += MUZ_ix_iz*( dzV * k_V_z[iz] + psi_V_z[idx] ) ;
  }else{ 
    Tz[in_idx] += MUZ_ix_iz*( dzV ) ;
  }

}
__global__ void cudaupdate_velocity(
	bool usetable,
	int * userestrict mat_index,
	bool usepunch,
	int current_step,
	int * userestrict first_arrive_step,
	int * userestrict last_step,
	float * userestrict Tx,
	float * userestrict Tz,
	float * userestrict V,

	float * userestrict  BV,

	float * userestrict psi_Tx_x,
	float * userestrict psi_Tz_z,

	float * userestrict b_Tx_x,
	float * userestrict b_Tz_z,

	float * userestrict c_Tx_x,
	float * userestrict c_Tz_z,

	float * userestrict k_Tx_x,
	float * userestrict k_Tz_z,
	
	int nx,
	int nz,
	int npml,
	int izshift)
{
  __shared__ float s_data[BDIMY+2*radius][BDIMX+2*radius];
  int ix = blockIdx.x*blockDim.x + threadIdx.x;
  int iz = blockIdx.y*blockDim.y + threadIdx.y+izshift;

  CHECK_IF_NEED_UPDATE;

  int ord;
  int idx;
  int in_idx = iz*nx + ix; // index for reading input
  int tx=threadIdx.x+radius;
  int ty=threadIdx.y+radius;
  float Tx_x,Tz_z;
  float BV_ix_iz;

  ZERO_HOLO;
  __syncthreads();

  if(usetable){
	BV_ix_iz= BV[mat_index[in_idx]];
  }else{
	BV_ix_iz= BV[in_idx];
  }

  /* stress imaging */
  if(iz==0) Tz[in_idx]=-Tz[in_idx+nx];


  ord=min(ix,iz);
//  ord=ix; // maybe top layer can also use 8th ,this reduce dispersion
  ord=min(ord,nx-1-ix);
  ord=min(ord,nz-1-iz);
  ord=min(ord,ORD8);

  EDGE_SHARE(Tx);
  __syncthreads();
  TxDIFF;
  __syncthreads();

  EDGE_SHARE(Tz);
  __syncthreads();
  TzDIFF;

  if(ix<npml){
    //left pml
    idx=ix+iz*2*npml;
    psi_Tx_x[idx] = b_Tx_x[ix] * psi_Tx_x[idx] + c_Tx_x[ix] * Tx_x;
    V[in_idx] += BV_ix_iz*( Tx_x * k_Tx_x[ix] + psi_Tx_x[idx] );
  }else if(ix>=nx-npml){
    //right pml
    idx=npml+nx-1-ix+iz*2*npml;
    psi_Tx_x[idx] = b_Tx_x[ix] * psi_Tx_x[idx] + c_Tx_x[ix] * Tx_x;
    V[in_idx] += BV_ix_iz*( Tx_x * k_Tx_x[ix] + psi_Tx_x[idx] );
  }else{
    V[in_idx] += BV_ix_iz*( Tx_x );
  }

  if(iz<npml && topabsorb){
	//top pml
	idx=(iz*nx)+ix;
    psi_Tz_z[idx] = b_Tz_z[iz] * psi_Tz_z[idx] + c_Tz_z[iz] * Tz_z;
    V[in_idx] += BV_ix_iz*(Tz_z * k_Tz_z[iz] + psi_Tz_z[idx] );

  }else if(iz>=nz-npml){
	// bottom
	idx=(npml+nz-1-iz)*nx+ix;
    psi_Tz_z[idx] = b_Tz_z[iz] * psi_Tz_z[idx] + c_Tz_z[iz] * Tz_z;
    V[in_idx] += BV_ix_iz*(Tz_z * k_Tz_z[iz] + psi_Tz_z[idx] );
  }else{
    V[in_idx] += BV_ix_iz*( Tz_z );
  }

}

void cu_step_stress(int deviceid,PARAM &param,FIELD & fld, MATERIAL &mat,CPML &cpml,SLIDE &slide)
{

   int nz1,nz2,tnz;
   int nx=fld.nx;
   int nz=fld.nz;
   int npml=cpml.npml;
   nz1=param.a_nz1[deviceid];
   nz2=param.a_nz2[deviceid];
   tnz=nz2-nz1+1;

   dim3 nblocks((nx+BDIMX-1)/BDIMX,(tnz+BDIMY-1)/BDIMY);
   dim3 blocksize(BDIMX,BDIMY);

   bool usepunch=slide.usepunch==1;
   if(mat.usetable){
   cudaupdate_stress<<<nblocks,blocksize>>>(
	mat.usetable,
	mat.index,
    usepunch, 
    slide.current_step,
    slide.first_arrive_step,
    slide.last_step,
	fld.Tx,
	fld.Tz,
	fld.V,
	mat.tbl_MUX,
	mat.tbl_MUZ,
	cpml.psi.V_x,
	cpml.psi.V_z,
	cpml.b.V_x,
    cpml.b.V_z,
	cpml.c.V_x,
    cpml.c.V_z,
	cpml.k.V_x,
    cpml.k.V_z,
	nx,
	nz,
	npml,
	nz1);
   }else{
   cudaupdate_stress<<<nblocks,blocksize>>>(
	mat.usetable,
	mat.index,
    usepunch, 
    slide.current_step,
    slide.first_arrive_step,
    slide.last_step,
	fld.Tx,
	fld.Tz,
	fld.V,
	mat.MUX,
	mat.MUZ,
	cpml.psi.V_x,
	cpml.psi.V_z,
	cpml.b.V_x,
    cpml.b.V_z,
	cpml.c.V_x,
    cpml.c.V_z,
	cpml.k.V_x,
    cpml.k.V_z,
	nx,
	nz,
	npml,
	nz1);
   }
   CUT_CHECK_ERROR("Error in step_forward_stress");
}

void cu_step_velocity(int deviceid,PARAM &param,FIELD & fld, MATERIAL &mat,CPML &cpml,SLIDE &slide)
{

   int nz1,nz2,tnz;
   int nx=fld.nx;
   int nz=fld.nz;
   int npml=cpml.npml;
   nz1=param.a_nz1[deviceid];
   nz2=param.a_nz2[deviceid];
   tnz=nz2-nz1+1;

   dim3 nblocks((nx+BDIMX-1)/BDIMX,(tnz+BDIMY-1)/BDIMY);
   dim3 blocksize(BDIMX,BDIMY);

   bool usepunch=slide.usepunch==1;
   if(mat.usetable){
   cudaupdate_velocity<<<nblocks,blocksize>>>(
	mat.usetable,
	mat.index,

    usepunch, 
    slide.current_step,
    slide.first_arrive_step,
    slide.last_step,

	fld.Tx,
	fld.Tz,
	fld.V,

	mat.tbl_BV,
	cpml.psi.Tx_x,
	cpml.psi.Tz_z,

	cpml.b.Tx_x,
	cpml.b.Tz_z,

	cpml.c.Tx_x,
	cpml.c.Tz_z,

	cpml.k.Tx_x,
	cpml.k.Tz_z,
	nx,
	nz,
	npml,
	nz1
	  );
   }else{
   cudaupdate_velocity<<<nblocks,blocksize>>>(
	mat.usetable,
	mat.index,

    usepunch, 
    slide.current_step,
    slide.first_arrive_step,
    slide.last_step,

	fld.Tx,
	fld.Tz,
	fld.V,

	mat.BV,

	cpml.psi.Tx_x,
	cpml.psi.Tz_z,

	cpml.b.Tx_x,
	cpml.b.Tz_z,

	cpml.c.Tx_x,
	cpml.c.Tz_z,

	cpml.k.Tx_x,
	cpml.k.Tz_z,
	nx,
	nz,
	npml,
	nz1
	  );
   };
   CUT_CHECK_ERROR("Error in step_forward_velocity");
}
