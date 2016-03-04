/* Author Dunzhu Li  dli@caltech.edu 
*/
#include"cpml.h"
#include"../gpu.h"
#include<cstdio>

CPML::CPML(int deviceid, CPML &cpml)
{

  npml=cpml.npml;
  nx=cpml.nx;
  nz=cpml.nz;
  pml_dt=cpml.pml_dt;
  pml_r=cpml.pml_r;
  pml_v=cpml.pml_v;
  pml_fc=cpml.pml_fc;

#define ALLOC_CPML(psi,comp,n)\
  safecall(cudaMalloc((void**)&psi.comp,sizeof(float)*n));\
  safecall(cudaMemcpy(psi.comp,cpml.psi.comp,sizeof(float)*n,cudaMemcpyHostToDevice));

  ALLOC_CPML(psi,Txx_x,2*npml*nz);
  ALLOC_CPML(psi,Txz_x,2*npml*nz);
  ALLOC_CPML(psi,  U_x,2*npml*nz);
  ALLOC_CPML(psi,  W_x,2*npml*nz);

  ALLOC_CPML(psi,Tzz_z,2*npml*nx);
  ALLOC_CPML(psi,Txz_z,2*npml*nx);
  ALLOC_CPML(psi,  U_z,2*npml*nx);
  ALLOC_CPML(psi,  W_z,2*npml*nx);

  ALLOC_CPML(b,Txx_x,nx);
  ALLOC_CPML(b,Txz_x,nx);
  ALLOC_CPML(b,  U_x,nx);
  ALLOC_CPML(b,  W_x,nx);

  ALLOC_CPML(b,Tzz_z,nz);
  ALLOC_CPML(b,Txz_z,nz);
  ALLOC_CPML(b,  U_z,nz);
  ALLOC_CPML(b,  W_z,nz);

  ALLOC_CPML(c,Txx_x,nx);
  ALLOC_CPML(c,Txz_x,nx);
  ALLOC_CPML(c,  U_x,nx);
  ALLOC_CPML(c,  W_x,nx);

  ALLOC_CPML(c,Tzz_z,nz);
  ALLOC_CPML(c,Txz_z,nz);
  ALLOC_CPML(c,  U_z,nz);
  ALLOC_CPML(c,  W_z,nz);

  ALLOC_CPML(k,Txx_x,nx);
  ALLOC_CPML(k,Txz_x,nx);
  ALLOC_CPML(k,  U_x,nx);
  ALLOC_CPML(k,  W_x,nx);

  ALLOC_CPML(k,Tzz_z,nz);
  ALLOC_CPML(k,Txz_z,nz);
  ALLOC_CPML(k,  U_z,nz);
  ALLOC_CPML(k,  W_z,nz);

  /*
   P is the position of \partial{U}/\partial{x}, U_x, 

   P-------U              P------U--------
           |                     |
           |                     |
           |                     |
           |                     |
      0:pml_pos              nx-1-pml_pos:nx-1

distance from boundary
      pml_pos-i               i-(nx-1-pml_pos)

   for example, nx=100 
   looks from U grid, pml boundary is at position 11,   0 : 11, 100-1-11 : 100-1
   looks from P grid, pml boundary is at position 11.5, 0 : 11, 

   W---------Txz
   |	      |
   |          |
   |          | 
   Txx(Tzz)---U ------ pml boundary
   			  |
			  |
			  |
			  |
			  pml boundary
   0------npml-0.5
   */

}

void CPML::cu_load_restart(int deviceid,CPML &cpml)
{
#define LOAD_CPML(psi,comp,n)\
  safecall(cudaMemcpy(psi.comp,cpml.psi.comp,sizeof(float)*n,cudaMemcpyHostToDevice));

  LOAD_CPML(psi,Txx_x,2*npml*nz);
  LOAD_CPML(psi,Txz_x,2*npml*nz);
  LOAD_CPML(psi,  U_x,2*npml*nz);
  LOAD_CPML(psi,  W_x,2*npml*nz);

  LOAD_CPML(psi,Tzz_z,2*npml*nx);
  LOAD_CPML(psi,Txz_z,2*npml*nx);
  LOAD_CPML(psi,  U_z,2*npml*nx);
  LOAD_CPML(psi,  W_z,2*npml*nx);

}

void CPML::cu_save_state(int deviceid,CPML &cpml)
{
#define CPYBACK_CPML(psi,comp,n)\
  safecall(cudaMemcpy(cpml.psi.comp,psi.comp,sizeof(float)*n,cudaMemcpyDeviceToHost));

  CPYBACK_CPML(psi,Txx_x,2*npml*nz);
  CPYBACK_CPML(psi,Txz_x,2*npml*nz);
  CPYBACK_CPML(psi,  U_x,2*npml*nz);
  CPYBACK_CPML(psi,  W_x,2*npml*nz);

  CPYBACK_CPML(psi,Tzz_z,2*npml*nx);
  CPYBACK_CPML(psi,Txz_z,2*npml*nx);
  CPYBACK_CPML(psi,  U_z,2*npml*nx);
  CPYBACK_CPML(psi,  W_z,2*npml*nx);

}
