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

#define s_data(i,j) s_data[j][i]

#define CHECK_IF_NEED_UPDATE  \
  if(usepunch){\
	if(\
		current_step <= first_arrive_step[(iz/BDIMY)*(nx/BDIMX)+(ix/BDIMX)] ||\
		current_step >= last_step        [(iz/BDIMY)*(nx/BDIMX)+(ix/BDIMX)] ){ \
	  return;\
	}\
  }


#define  UDIFF \
  ux=d_coef[ord][0]*( s_data(tx,ty)   - s_data(tx-1,ty) )\
    -d_coef[ord][1]*( s_data(tx+1,ty) - s_data(tx-2,ty) )\
    +d_coef[ord][2]*( s_data(tx+2,ty) - s_data(tx-3,ty) )\
    -d_coef[ord][3]*( s_data(tx+3,ty) - s_data(tx-4,ty) );\
  uz=d_coef[ord][0]*( s_data(tx,ty)   - s_data(tx,ty-1) )\
    -d_coef[ord][1]*( s_data(tx,ty+1) - s_data(tx,ty-2) )\
    +d_coef[ord][2]*( s_data(tx,ty+2) - s_data(tx,ty-3) )\
    -d_coef[ord][3]*( s_data(tx,ty+3) - s_data(tx,ty-4) );

#define WDIFF \
  wx=d_coef[ord][0]*( s_data(tx+1,ty) - s_data(tx,ty) ) \
    -d_coef[ord][1]*( s_data(tx+2,ty) - s_data(tx-1,ty) ) \
    +d_coef[ord][2]*( s_data(tx+3,ty) - s_data(tx-2,ty) ) \
    -d_coef[ord][3]*( s_data(tx+4,ty) - s_data(tx-3,ty) );\
  wz=d_coef[ord][0]*( s_data(tx,ty+1) - s_data(tx,ty) )\
    -d_coef[ord][1]*( s_data(tx,ty+2) - s_data(tx,ty-1) )\
    +d_coef[ord][2]*( s_data(tx,ty+3) - s_data(tx,ty-2) )\
    -d_coef[ord][3]*( s_data(tx,ty+4) - s_data(tx,ty-3) );

#define  DXTXX \
  dxTxx=d_coef[ord][0]*( s_data(tx+1,ty)   - s_data(tx,ty) )\
       -d_coef[ord][1]*( s_data(tx+2,ty) - s_data(tx-1,ty) )\
       +d_coef[ord][2]*( s_data(tx+3,ty) - s_data(tx-2,ty) )\
       -d_coef[ord][3]*( s_data(tx+4,ty) - s_data(tx-3,ty) );


#define  DZTZZ \
  dzTzz=d_coef[ord][0]*( s_data(tx,ty)   - s_data(tx,ty-1) )\
       -d_coef[ord][1]*( s_data(tx,ty+1) - s_data(tx,ty-2) )\
       +d_coef[ord][2]*( s_data(tx,ty+2) - s_data(tx,ty-3) )\
       -d_coef[ord][3]*( s_data(tx,ty+3) - s_data(tx,ty-4) );


#define DXDZTXZ \
  dxTxz=d_coef[ord][0]*( s_data(tx,ty)   - s_data(tx-1,ty) ) \
       -d_coef[ord][1]*( s_data(tx+1,ty) - s_data(tx-2,ty) ) \
       +d_coef[ord][2]*( s_data(tx+2,ty) - s_data(tx-3,ty) ) \
       -d_coef[ord][3]*( s_data(tx+3,ty) - s_data(tx-4,ty) );\
  dzTxz=d_coef[ord][0]*( s_data(tx,ty+1) - s_data(tx,ty) )\
       -d_coef[ord][1]*( s_data(tx,ty+2) - s_data(tx,ty-1) )\
       +d_coef[ord][2]*( s_data(tx,ty+3) - s_data(tx,ty-2) )\
       -d_coef[ord][3]*( s_data(tx,ty+4) - s_data(tx,ty-3) );

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
	int * __restrict__ mat_index,
	bool usepunch,
	int current_step,
	int * __restrict__ first_arrive_step,
	int * __restrict__ last_step,
	float * __restrict__ Txx,
	float * __restrict__ Txz,
	float * __restrict__ Tzz,
	float * __restrict__ U,
	float * __restrict__ W,
	float * __restrict__  MU,
	float * __restrict__ MUA,
	float * __restrict__ LAM,
	float * __restrict__ psi_U_x,
	float * __restrict__ psi_W_x,
	float * __restrict__ psi_U_z,
	float * __restrict__ psi_W_z,
	float * __restrict__ b_U_x,
    float * __restrict__ b_W_x,
    float * __restrict__ b_U_z,
    float * __restrict__ b_W_z,
	float * __restrict__ c_U_x,
    float * __restrict__ c_W_x,
    float * __restrict__ c_U_z,
    float * __restrict__ c_W_z,
	float * __restrict__ k_U_x,
    float * __restrict__ k_W_x,
    float * __restrict__ k_U_z,
    float * __restrict__ k_W_z,
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
  float ux,uz,wx,wz;
  float LAM_ix_iz;
  float MU_ix_iz ;
  float MUA_ix_iz;
  float gam     ; 

  if(usetable){
	LAM_ix_iz= LAM[mat_index[in_idx]];
	MU_ix_iz = MU [mat_index[in_idx]];
	MUA_ix_iz= MUA[mat_index[in_idx]];
  }else{
	LAM_ix_iz= LAM[in_idx];
	MU_ix_iz = MU [in_idx];
	MUA_ix_iz= MUA[in_idx];
  }

  gam      = LAM_ix_iz + 2.0f*MU_ix_iz;

  ZERO_HOLO;
  __syncthreads();

    
  ord=min(ix,iz);
//  ord=ix; // maybe top layer can also use 8th,this reduced dispersion
  ord=min(ord,nx-1-ix);
  ord=min(ord,nz-1-iz);
  ord=min(ord,ORD8);

  EDGE_SHARE(U);
  __syncthreads();
  UDIFF; 
  __syncthreads();

  EDGE_SHARE(W);
  __syncthreads();
  WDIFF;
#ifdef NOPML
  if(iz==0){ // free edge
	Txx[in_idx]+= 2.0f*MU_ix_iz*ux;
	Tzz[in_idx]= 0.0;
	Txz[in_idx]= 0.0;
  }else{
	Txx[in_idx] += ( gam*ux + LAM_ix_iz*wz );
	Tzz[in_idx] += ( LAM_ix_iz*ux + gam*wz );
	Txz[in_idx] += ( MUA_ix_iz*(uz + wx)   );
  }
#else
  if(ix<npml){
    //left pml
    idx=ix+iz*2*npml;
    psi_U_x[idx] = b_U_x[ix] * psi_U_x[idx] + c_U_x[ix] * ux;
    Txx[in_idx] += gam*( ux * k_U_x[ix] + psi_U_x[idx] );
    Tzz[in_idx] += LAM_ix_iz*( ux * k_U_x[ix] + psi_U_x[idx] );

    psi_W_x[idx] = b_W_x[ix] * psi_W_x[idx] + c_W_x[ix] * wx;
    Txz[in_idx] += MUA_ix_iz*( wx * k_W_x[ix] + psi_W_x[idx] );

  }else if(ix>=nx-npml){
    //right pml
    idx=npml+nx-1-ix+iz*2*npml;
    psi_U_x[idx] = b_U_x[ix] * psi_U_x[idx] + c_U_x[ix] * ux;
    Txx[in_idx] += gam*( ux * k_U_x[ix] + psi_U_x[idx] );
    Tzz[in_idx] += LAM_ix_iz*( ux * k_U_x[ix] + psi_U_x[idx] );

    psi_W_x[idx] = b_W_x[ix] * psi_W_x[idx] + c_W_x[ix] * wx;
    Txz[in_idx] += MUA_ix_iz*( wx * k_W_x[ix] + psi_W_x[idx] );
  }else if(iz!=0){
    Txx[in_idx] += gam*( ux );
    Tzz[in_idx] += LAM_ix_iz*( ux );
    Txz[in_idx] += MUA_ix_iz*( wx );
  }

  if(iz<npml && topabsorb){
	//top pml
	idx=(iz*nx)+ix;
    psi_W_z[idx] = b_W_z[iz] * psi_W_z[idx] + c_W_z[iz] * wz;
    Txx[in_idx] += LAM_ix_iz*( wz * k_W_z[iz] + psi_W_z[idx] ) ;
    Tzz[in_idx] += gam      *( wz * k_W_z[iz] + psi_W_z[idx] ) ;

    psi_U_z[idx] = b_U_z[iz] * psi_U_z[idx] + c_U_z[iz] * uz;
    Txz[in_idx] += MUA_ix_iz*( uz * k_U_z[iz] + psi_U_z[idx] ) ;
  }else if(iz>=nz-npml){
	// bottom
	idx=(npml+nz-1-iz)*nx+ix;
    psi_W_z[idx] = b_W_z[iz] * psi_W_z[idx] + c_W_z[iz] * wz;
    Txx[in_idx] += LAM_ix_iz*( wz * k_W_z[iz] + psi_W_z[idx] ) ;
    Tzz[in_idx] += gam      *( wz * k_W_z[iz] + psi_W_z[idx] ) ;

    psi_U_z[idx] = b_U_z[iz] * psi_U_z[idx] + c_U_z[iz] * uz;
    Txz[in_idx] += MUA_ix_iz*( uz * k_U_z[iz] + psi_U_z[idx] ) ;
  }else if(iz==0){ 
	// lam=0,gam=2*mu, so gam*ux+lam*wz=2*mu*ux
	Txx[in_idx]+= 2.0f*MU_ix_iz*ux;
	Tzz[in_idx]= 0.0;
	Txz[in_idx]= 0.0;
  }else{
    Txx[in_idx] += LAM_ix_iz*( wz );
    Tzz[in_idx] += gam      *( wz );
    Txz[in_idx] += MUA_ix_iz*( uz );
  }
#endif

}
__global__ void cudaupdate_velocity(
	bool usetable,
	int * __restrict__ mat_index,
	bool usepunch,
	int current_step,
	int * __restrict__ first_arrive_step,
	int * __restrict__ last_step,
	float * __restrict__ Txx,
	float * __restrict__ Txz,
	float * __restrict__ Tzz,
	float * __restrict__ U,
	float * __restrict__ W,

	float * __restrict__  BU,
	float * __restrict__  BW,

	float * __restrict__ psi_Txx_x,
	float * __restrict__ psi_Txz_x,
	float * __restrict__ psi_Txz_z,
	float * __restrict__ psi_Tzz_z,

	float * __restrict__ b_Txx_x,
	float * __restrict__ b_Txz_x,
	float * __restrict__ b_Txz_z,
	float * __restrict__ b_Tzz_z,

	float * __restrict__ c_Txx_x,
	float * __restrict__ c_Txz_x,
	float * __restrict__ c_Txz_z,
	float * __restrict__ c_Tzz_z,

	float * __restrict__ k_Txx_x,
	float * __restrict__ k_Txz_x,
	float * __restrict__ k_Txz_z,
	float * __restrict__ k_Tzz_z,
	
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
  float dxTxx,dzTzz,dxTxz,dzTxz;
  float BU_ix_iz;
  float BW_ix_iz;

  ZERO_HOLO;
  __syncthreads();

  if(usetable){
	BU_ix_iz= BU[mat_index[in_idx]];
	BW_ix_iz= BW[mat_index[in_idx]];
  }else{
	BU_ix_iz= BU[in_idx];
	BW_ix_iz= BW[in_idx];
  }


  ord=min(ix,iz);
//  ord=ix; // maybe top layer can also use 8th ,this reduce dispersion
  ord=min(ord,nx-1-ix);
  ord=min(ord,nz-1-iz);
  ord=min(ord,ORD8);

  EDGE_SHARE(Txx);
  __syncthreads();
  DXTXX;
  __syncthreads();

  EDGE_SHARE(Tzz);
  __syncthreads();
  DZTZZ;
  __syncthreads();

  EDGE_SHARE(Txz);
  __syncthreads();
  DXDZTXZ;
#ifdef NOPML
  if(iz==0){
	U[in_idx] += BU_ix_iz*(dxTxx+dzTxz);
	W[in_idx]=0.0;
  }else{
    U[in_idx] += BU_ix_iz*(dxTxx+dzTxz);
    W[in_idx] += BW_ix_iz*(dxTxz+dzTzz);
  }
#else
  if(ix<npml){
    //left pml
    idx=ix+iz*2*npml;
    psi_Txx_x[idx] = b_Txx_x[ix] * psi_Txx_x[idx] + c_Txx_x[ix] * dxTxx;
    U[in_idx] += BU_ix_iz*( dxTxx * k_Txx_x[ix] + psi_Txx_x[idx] );

    psi_Txz_x[idx] = b_Txz_x[ix] * psi_Txz_x[idx] + c_Txz_x[ix] * dxTxz;
    W[in_idx] += BW_ix_iz*( dxTxz * k_Txz_x[ix] + psi_Txz_x[idx] );
  }else if(ix>=nx-npml){
    //right pml
    idx=npml+nx-1-ix+iz*2*npml;
    psi_Txx_x[idx] = b_Txx_x[ix] * psi_Txx_x[idx] + c_Txx_x[ix] * dxTxx;
    U[in_idx] += BU_ix_iz*( dxTxx * k_Txx_x[ix] + psi_Txx_x[idx] );

    psi_Txz_x[idx] = b_Txz_x[ix] * psi_Txz_x[idx] + c_Txz_x[ix] * dxTxz;
    W[in_idx] += BW_ix_iz*( dxTxz * k_Txz_x[ix] + psi_Txz_x[idx] );
  }else if(iz!=0){
    U[in_idx] += BU_ix_iz*(dxTxx);
    W[in_idx] += BW_ix_iz*(dxTxz);
  }

  if(iz<npml && topabsorb){
	//top pml
	idx=(iz*nx)+ix;
    psi_Txz_z[idx] = b_Txz_z[iz] * psi_Txz_z[idx] + c_Txz_z[iz] * dzTxz;
    psi_Tzz_z[idx] = b_Tzz_z[iz] * psi_Tzz_z[idx] + c_Tzz_z[iz] * dzTzz;

    U[in_idx] += BU_ix_iz*(dzTxz * k_Txz_z[iz] + psi_Txz_z[idx] );
    W[in_idx] += BW_ix_iz*(dzTzz * k_Tzz_z[iz] + psi_Tzz_z[idx] );

  }else if(iz>=nz-npml){
	// bottom
	idx=(npml+nz-1-iz)*nx+ix;
    psi_Txz_z[idx] = b_Txz_z[iz] * psi_Txz_z[idx] + c_Txz_z[iz] * dzTxz;
    psi_Tzz_z[idx] = b_Tzz_z[iz] * psi_Tzz_z[idx] + c_Tzz_z[iz] * dzTzz;

    U[in_idx] += BU_ix_iz*(dzTxz * k_Txz_z[iz] + psi_Txz_z[idx] );
    W[in_idx] += BW_ix_iz*(dzTzz * k_Tzz_z[iz] + psi_Tzz_z[idx] );
  }else if(iz==0){
	U[in_idx] += BU_ix_iz*(dxTxx+dzTxz);
	W[in_idx]=0.0;
  }else{
    U[in_idx] += BU_ix_iz*(dzTxz);
    W[in_idx] += BW_ix_iz*(dzTzz);
  }
#endif
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
	fld.Txx,
	fld.Txz,
	fld.Tzz,
	fld.U,
	fld.W,
	mat.tbl_MU,
	mat.tbl_MUA,
	mat.tbl_LAM,
	cpml.psi.U_x,
	cpml.psi.W_x,
	cpml.psi.U_z,
	cpml.psi.W_z,
	cpml.b.U_x,
    cpml.b.W_x,
    cpml.b.U_z,
    cpml.b.W_z,
	cpml.c.U_x,
    cpml.c.W_x,
    cpml.c.U_z,
    cpml.c.W_z,
	cpml.k.U_x,
    cpml.k.W_x,
    cpml.k.U_z,
    cpml.k.W_z,
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

	fld.Txx,
	fld.Txz,
	fld.Tzz,
	fld.U,
	fld.W,

	mat.MU,
	mat.MUA,
	mat.LAM,
	cpml.psi.U_x,
	cpml.psi.W_x,
	cpml.psi.U_z,
	cpml.psi.W_z,
	cpml.b.U_x,
    cpml.b.W_x,
    cpml.b.U_z,
    cpml.b.W_z,
	cpml.c.U_x,
    cpml.c.W_x,
    cpml.c.U_z,
    cpml.c.W_z,
	cpml.k.U_x,
    cpml.k.W_x,
    cpml.k.U_z,
    cpml.k.W_z,
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

	fld.Txx,
	fld.Txz,
	fld.Tzz,
	fld.U,
	fld.W,

	mat.tbl_BU,
	mat.tbl_BW,

	cpml.psi.Txx_x,
	cpml.psi.Txz_x,
	cpml.psi.Txz_z,
	cpml.psi.Tzz_z,

	cpml.b.Txx_x,
	cpml.b.Txz_x,
	cpml.b.Txz_z,
	cpml.b.Tzz_z,

	cpml.c.Txx_x,
	cpml.c.Txz_x,
	cpml.c.Txz_z,
	cpml.c.Tzz_z,

	cpml.k.Txx_x,
	cpml.k.Txz_x,
	cpml.k.Txz_z,
	cpml.k.Tzz_z,
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

	fld.Txx,
	fld.Txz,
	fld.Tzz,
	fld.U,
	fld.W,

	mat.BU,
	mat.BW,

	cpml.psi.Txx_x,
	cpml.psi.Txz_x,
	cpml.psi.Txz_z,
	cpml.psi.Tzz_z,

	cpml.b.Txx_x,
	cpml.b.Txz_x,
	cpml.b.Txz_z,
	cpml.b.Tzz_z,

	cpml.c.Txx_x,
	cpml.c.Txz_x,
	cpml.c.Txz_z,
	cpml.c.Tzz_z,

	cpml.k.Txx_x,
	cpml.k.Txz_x,
	cpml.k.Txz_z,
	cpml.k.Tzz_z,
	nx,
	nz,
	npml,
	nz1
	  );
   };
   CUT_CHECK_ERROR("Error in step_forward_velocity");
}