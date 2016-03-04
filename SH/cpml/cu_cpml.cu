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

  ALLOC_CPML(psi,Tx_x,2*npml*nz);
  ALLOC_CPML(psi, V_x,2*npml*nz);

  ALLOC_CPML(psi,Tz_z,2*npml*nx);
  ALLOC_CPML(psi, V_z,2*npml*nx);

  ALLOC_CPML(b,Tx_x,nx);
  ALLOC_CPML(b, V_x,nx);

  ALLOC_CPML(b,Tz_z,nz);
  ALLOC_CPML(b, V_z,nz);

  ALLOC_CPML(c,Tx_x,nx);
  ALLOC_CPML(c, V_x,nx);

  ALLOC_CPML(c,Tz_z,nz);
  ALLOC_CPML(c, V_z,nz);

  ALLOC_CPML(k,Tx_x,nx);
  ALLOC_CPML(k, V_x,nx);

  ALLOC_CPML(k,Tz_z,nz);
  ALLOC_CPML(k, V_z,nz);

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
