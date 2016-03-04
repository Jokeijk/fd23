#include"cpml.h"
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include"../param/param.h"

int helper_pmlparameter(int npml,float pml_pos,float startpos,int nx,float *b,float*c,float*k,float dt,float R,float v,float fc)
{
  float *sg_u =new float[nx];
  float *kp_u =new float[nx];
  float *al_u =new float[nx];
  float *b_u  =new float[nx];
  float *c_u  =new float[nx];
  float *d    =new float[nx];

  const float pk=2.0;
  const float psig=2.0;
  const float palpha=1.0;

  float sigmax=-(psig+1.0)/(2.0*floor(pml_pos))*v*log(R);
  float alphamax=3.141592653589*fc;
  float kmax=2.0;

  for(int i=0;i<nx;i++){
	sg_u [i]=0.0;
	kp_u [i]=1.0;
	al_u [i]=alphamax;
	d    [i]=0.0;
  }
  for(int i=0;i<npml+2;i++){
	float tmp=( pml_pos-(i+startpos) )/floor(npml);
	if(tmp>0)	  d[i]=tmp;
  }
  for(int i=nx-npml-2;i<nx;i++){
	float tmp=( (i+startpos)-(nx-1-pml_pos) )/floor(npml);
	if(tmp>0)	  d[i]=tmp;
  }

  for(int i=0;i<nx;i++){
	sg_u[i]=sigmax * pow(d[i], psig);
	kp_u[i]=1.0+(kmax-1.0)*pow(d[i], pk);
	al_u[i]=alphamax *(1.0- pow(d[i],palpha));
  }

  for(int i=0;i<nx;i++){
	b_u [i]=exp(-(sg_u[i]/kp_u[i] +  al_u[i])*dt);
	if(sg_u[i]+kp_u[i]*al_u[i]<1E-15){
	  c_u[i]=0;
	}else{
	  c_u [i]=sg_u[i]*(b_u[i]-1.0)/(sg_u[i]+kp_u[i]*al_u[i])/kp_u[i];
	}
	kp_u[i]=1.0/kp_u[i];
  }
  //debug, forget why I did this
  for(int i=0;i<4;i++){
	b_u[i]=0.0;
	c_u[i]=0.0;
	kp_u[i]=0.0;
  }
  for(int i=nx-5;i<nx;i++){
	b_u[i]=0.0;
	c_u[i]=0.0;
	kp_u[i]=0.0;
  }
  //end


  std::copy(b_u,b_u+nx,b);
  std::copy(c_u,c_u+nx,c);
  std::copy(kp_u,kp_u+nx,k);

//  safecall(cudaMemcpy(b, b_u,sizeof(float)*nx,cudaMemcpyHostToDevice));
//  safecall(cudaMemcpy(c, c_u,sizeof(float)*nx,cudaMemcpyHostToDevice));
//  safecall(cudaMemcpy(k,kp_u,sizeof(float)*nx,cudaMemcpyHostToDevice));

  delete []sg_u; 
  delete []kp_u; 
  delete []al_u; 
  delete []b_u ; 
  delete []c_u ; 
  delete []d   ; 
  return 0;
}
CPML::CPML(PARAM &param)
{
  int npml=param.npml;
  int nx=param.nx;
  int nz=param.nz;
  float pml_dt=param.pml_dt;
  float pml_r=param.pml_r;
  float pml_v=param.pml_v;
  float pml_fc=param.pml_fc;

  this->npml=npml;
  this->nx=nx;
  this->nz=nz;
  this->pml_dt=pml_dt;
  this->pml_r=pml_r;
  this->pml_v=pml_v;
  this->pml_fc=pml_fc;

  float dt=pml_dt;

#define ALLOC_CPML(psi,comp,n)\
  psi.comp=new float[n];std::fill_n(psi.comp,n,0);

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

  helper_pmlparameter(npml,npml-0.5 ,0.5,nx,b.W_x  ,c.W_x  ,k.W_x  ,pml_dt,pml_r,pml_v,pml_fc);
  helper_pmlparameter(npml,npml-0.5 ,0.5,nx,b.Txx_x,c.Txx_x,k.Txx_x,pml_dt,pml_r,pml_v,pml_fc);
  helper_pmlparameter(npml,npml-0.5 ,0.0,nx,b.U_x  ,c.U_x  ,k.U_x  ,pml_dt,pml_r,pml_v,pml_fc);
  helper_pmlparameter(npml,npml-0.5 ,0.0,nx,b.Txz_x,c.Txz_x,k.Txz_x,pml_dt,pml_r,pml_v,pml_fc);

  helper_pmlparameter(npml,npml-0.5 ,0.5,nz,b.W_z  ,c.W_z  ,k.W_z  ,pml_dt,pml_r,pml_v,pml_fc);
  helper_pmlparameter(npml,npml-0.5 ,0.5,nz,b.Txz_z,c.Txz_z,k.Txz_z,pml_dt,pml_r,pml_v,pml_fc);
  helper_pmlparameter(npml,npml-0.5 ,0.0,nz,b.U_z  ,c.U_z  ,k.U_z  ,pml_dt,pml_r,pml_v,pml_fc);
  helper_pmlparameter(npml,npml-0.5 ,0.0,nz,b.Tzz_z,c.Tzz_z,k.Tzz_z,pml_dt,pml_r,pml_v,pml_fc);
}
