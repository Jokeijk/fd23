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
  //debug
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

   Tz---------
   |	      |
   |          |
   |          | 
   V---------Tx ------ pml boundary
   			  |
			  |
			  |
			  |
			  pml boundary
   0------npml-0.5
   */
  helper_pmlparameter(npml,npml-0.5 ,0.0,nx,b.Tx_x, c.Tx_x, k.Tx_x  ,pml_dt,pml_r,pml_v,pml_fc);
  helper_pmlparameter(npml,npml-0.5 ,0.5,nz,b.Tz_z, c.Tz_z, k.Tz_z  ,pml_dt,pml_r,pml_v,pml_fc);

  helper_pmlparameter(npml,npml-0.5 ,0.5,nx,b.V_x, c.V_x, k.V_x,pml_dt,pml_r,pml_v,pml_fc);
  helper_pmlparameter(npml,npml-0.5 ,0.0,nz,b.V_z, c.V_z, k.V_z,pml_dt,pml_r,pml_v,pml_fc);

}
